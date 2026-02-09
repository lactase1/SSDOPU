#!/usr/bin/env python3
"""
dcm_to_mp4.py

Convert DICOM (.dcm) multi-frame files to MP4 video(s).

Usage:
    python dcm_to_mp4.py <input_path> [--output <out_path_or_dir>] [--fps FPS] [--verbose]

If <input_path> is a file, produces one MP4. If it's a directory, processes all .dcm files inside and
writes per-file MP4s to the output directory (or same folder if not given).

Dependencies: pydicom, opencv-python, numpy

Examples:
    python dcm_to_mp4.py data/example.dcm --output out.mp4 --fps 10
    python dcm_to_mp4.py data/dcm_folder --output out_dir --fps 15
"""

import os
import argparse
import numpy as np
import cv2
import pydicom

# Try to import imageio for better MP4 support
try:
    import imageio
    HAS_IMAGEIO = True
except ImportError:
    HAS_IMAGEIO = False


def get_frames_from_dataset(ds, verbose=False):
    """Return list of frames as numpy arrays with shape (rows, cols, [channels])"""
    # Try pydicom's pixel_array first
    try:
        arr = ds.pixel_array
        if verbose:
            print('Loaded pixel_array, shape=', arr.shape, 'dtype=', arr.dtype)
    except Exception as e:
        if verbose:
            print('ds.pixel_array failed, trying manual PixelData handling:', e)
        # Manual fallback: use PixelData unpacking for bit-packed formats
        rows = int(getattr(ds, 'Rows', 1))
        cols = int(getattr(ds, 'Columns', 1))
        frames = int(getattr(ds, 'NumberOfFrames', 1))
        samples = int(getattr(ds, 'SamplesPerPixel', 1))
        bits = int(getattr(ds, 'BitsAllocated', 8))
        raw = ds.PixelData
        if bits == 1:
            data_bytes = np.frombuffer(raw, dtype=np.uint8)
            all_bits = np.unpackbits(data_bytes)
            needed = rows * cols * frames * samples
            if all_bits.size < needed:
                raise ValueError(f'PixelData bits not enough: need {needed}, got {all_bits.size}')
            all_bits = all_bits[:needed]
            if samples == 1:
                arr = all_bits.reshape((frames, rows, cols)).transpose(1,2,0)
            else:
                arr = all_bits.reshape((frames, rows, cols, samples)).transpose(1,2,0,3)
        else:
            if bits <= 8:
                dtype = np.uint8
            elif bits <= 16:
                dtype = np.uint16
            else:
                dtype = np.uint32
            arr_flat = np.frombuffer(raw, dtype=dtype)
            expected = rows * cols * frames * samples
            if arr_flat.size < expected:
                raise ValueError(f'PixelData length not enough: need {expected}, got {arr_flat.size}')
            arr_flat = arr_flat[:expected]
            if frames > 1:
                if samples == 1:
                    arr = arr_flat.reshape((frames, rows, cols)).transpose(1,2,0)
                else:
                    arr = arr_flat.reshape((frames, rows, cols, samples)).transpose(1,2,0,3)
            else:
                if samples == 1:
                    arr = arr_flat.reshape((rows, cols))
                else:
                    arr = arr_flat.reshape((rows, cols, samples))
        if verbose:
            print('Manual PixelData parsed, arr.shape=', arr.shape, 'dtype=', arr.dtype)

    # Normalize arr into frame list
    frames_list = []
    # Determine number of frames
    if hasattr(ds, 'NumberOfFrames') and int(getattr(ds, 'NumberOfFrames')) > 1:
        nframes = int(getattr(ds, 'NumberOfFrames'))
    else:
        # try to infer
        nframes = 1
        for dim in arr.shape:
            if dim > 1 and dim != getattr(ds, 'Rows', dim) and dim != getattr(ds, 'Columns', dim):
                nframes = max(nframes, dim)
        # fallback: if arr has 3 dims and last dim not 3 (not RGB), assume it's frames
        if arr.ndim == 3 and arr.shape[2] != 3:
            nframes = arr.shape[2]

    # Move frames axis to front if needed
    if arr.ndim == 2:
        frames_list = [arr]
    elif arr.ndim == 3:
        # Could be (rows, cols, frames) or (frames, rows, cols) or (rows, cols, channels)
        if 'PhotometricInterpretation' in ds and 'RGB' in ds.PhotometricInterpretation and arr.shape[2] == 3:
            # RGB single frame
            frames_list = [arr]
        elif arr.shape[0] == nframes:
            # (frames, rows, cols)
            for i in range(arr.shape[0]):
                frames_list.append(arr[i])
        elif arr.shape[2] == nframes:
            # (rows, cols, frames)
            for i in range(arr.shape[2]):
                frames_list.append(arr[:, :, i])
        else:
            # assume frames along axis 2 if it's large
            for i in range(arr.shape[2]):
                frames_list.append(arr[:, :, i])
    elif arr.ndim == 4:
        # Most likely (rows, cols, samples, frames) or (frames, rows, cols, samples)
        if arr.shape[3] == nframes:
            # (rows, cols, samples, frames)
            for i in range(arr.shape[3]):
                frame = arr[:, :, :, i]
                frames_list.append(frame)
        elif arr.shape[0] == nframes:
            # (frames, rows, cols, samples)
            for i in range(arr.shape[0]):
                frames_list.append(arr[i])
        else:
            # fallback: iterate last axis
            for i in range(arr.shape[-1]):
                frames_list.append(arr.take(i, axis=-1))
    else:
        raise ValueError('Unsupported pixel array dimensionality: %s' % (arr.shape,))

    return frames_list


def normalize_frame_to_uint8(frame):
    """Normalize numpy frame to uint8 2D or 3D array suitable for video writing"""
    if frame.dtype == np.uint8:
        return frame
    if frame.dtype == np.bool_:
        return (frame.astype(np.uint8) * 255)
    fd = frame.astype(np.float64)
    mn = fd.min()
    mx = fd.max()
    if mx == mn:
        # constant image -> expand to 64x64
        val = 0
        if not np.isnan(mn):
            val = int(round((mn - mn) / max(1, mx - mn) * 255))
        return np.full((64, 64), val, dtype=np.uint8)
    norm = (fd - mn) / (mx - mn)
    scaled = (norm * 255.0).clip(0, 255).astype(np.uint8)
    return scaled


def frames_to_mp4(frames, out_path, fps=10, verbose=False, skip_existing=False):
    """Save frames as MP4 video using imageio (preferred) or OpenCV fallback"""
    if len(frames) == 0:
        raise ValueError('No frames to save')
    
    # Prepare first frame to get dimensions
    first_frame = normalize_frame_to_uint8(frames[0])
    if first_frame.ndim == 2:
        height, width = first_frame.shape
    else:
        height, width = first_frame.shape[:2]
    
    # Ensure output path has .mp4 extension
    out_path = os.path.splitext(out_path)[0] + '.mp4'
    
    # Prepare all frames as RGB uint8
    rgb_frames = []
    for f in frames:
        f_u8 = normalize_frame_to_uint8(f)
        if f_u8.ndim == 2:
            rgb_frame = cv2.cvtColor(f_u8, cv2.COLOR_GRAY2RGB)
        elif f_u8.ndim == 3 and f_u8.shape[2] == 3:
            rgb_frame = f_u8  # Already RGB
        elif f_u8.ndim == 3 and f_u8.shape[2] == 4:
            rgb_frame = cv2.cvtColor(f_u8, cv2.COLOR_RGBA2RGB)
        else:
            rgb_frame = cv2.cvtColor(f_u8[..., 0], cv2.COLOR_GRAY2RGB)
        rgb_frames.append(rgb_frame)
    
    # Try imageio first (more reliable for MP4)
    if HAS_IMAGEIO:
        try:
            if verbose:
                print(f'Using imageio to write MP4: {out_path}')
            
            # Check if file exists and try to remove it
            if os.path.exists(out_path):
                try:
                    os.remove(out_path)
                except PermissionError:
                    if skip_existing:
                        if verbose:
                            print(f'Skipping {out_path} - file exists and may be in use')
                        return None
                    raise PermissionError(f'Cannot overwrite {out_path} - file may be in use')
            
            # Check write permission by testing parent directory
            out_dir = os.path.dirname(out_path) or '.'
            if not os.access(out_dir, os.W_OK):
                raise PermissionError(f'No write permission for directory: {out_dir}')
            
            # Use imageio with ffmpeg backend
            writer = imageio.get_writer(out_path, fps=fps, codec='libx264', 
                                        pixelformat='yuv420p', macro_block_size=1)
            for frame in rgb_frames:
                writer.append_data(frame)
            writer.close()
            if verbose:
                print(f'Video saved: {out_path} ({len(frames)} frames, fps={fps})')
            return out_path
        except PermissionError as e:
            raise e  # Re-raise permission errors
        except Exception as e:
            if verbose:
                print(f'imageio failed: {e}, trying OpenCV fallback...')
    
    # OpenCV fallback with multiple codecs
    codecs = [
        ('avc1', '.mp4'),
        ('mp4v', '.mp4'),
        ('X264', '.mp4'),
        ('H264', '.mp4'),
    ]
    
    writer = None
    final_out_path = out_path
    for codec, ext in codecs:
        final_out_path = os.path.splitext(out_path)[0] + ext
        
        fourcc = cv2.VideoWriter_fourcc(*codec)
        writer = cv2.VideoWriter(final_out_path, fourcc, fps, (width, height))
        
        if writer.isOpened():
            if verbose:
                print(f'Using OpenCV codec: {codec}, output: {final_out_path}')
            break
        else:
            writer.release()
            writer = None
            if verbose:
                print(f'Codec {codec} failed, trying next...')
    
    if writer is None or not writer.isOpened():
        raise ValueError(f'Failed to open video writer for {out_path}. No suitable codec found. '
                        f'Try: pip install imageio imageio-ffmpeg')
    
    for i, f in enumerate(rgb_frames):
        # Convert RGB to BGR for OpenCV
        bgr_frame = cv2.cvtColor(f, cv2.COLOR_RGB2BGR)
        writer.write(bgr_frame)
        if verbose:
            print(f'Written frame {i+1}/{len(frames)} size=({width}, {height})')
    
    writer.release()
    if verbose:
        print(f'Video saved: {final_out_path} ({len(frames)} frames, fps={fps})')
    
    return final_out_path


def process_file(dcm_path, out_path=None, fps=10, verbose=False, skip_existing=False):
    ds = pydicom.dcmread(dcm_path)
    frames = get_frames_from_dataset(ds, verbose=verbose)
    if out_path is None:
        out_path = os.path.splitext(dcm_path)[0] + '.mp4'
    elif os.path.isdir(out_path):
        out_path = os.path.join(out_path, os.path.splitext(os.path.basename(dcm_path))[0] + '.mp4')

    final_out_path = frames_to_mp4(frames, out_path, fps=fps, verbose=verbose, skip_existing=skip_existing)
    return final_out_path


def main():
    parser = argparse.ArgumentParser(description='Convert DICOM (multi-frame) to MP4 video')
    parser.add_argument('input', help='DICOM file or directory containing .dcm files')
    parser.add_argument('--output', '-o', help='Output MP4 file path or directory (defaults to input folder)', default=None)
    parser.add_argument('--fps', '-f', type=int, default=10, help='Frames per second (default: 10)')
    parser.add_argument('--skip-existing', '-s', action='store_true', help='Skip files that already exist and cannot be overwritten')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    if os.path.isdir(args.input):
        out_dir = args.output if args.output is not None else args.input
        os.makedirs(out_dir, exist_ok=True)
        dcm_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.lower().endswith('.dcm')]
        if not dcm_files:
            print('No .dcm files found in', args.input)
            return
        for f in dcm_files:
            try:
                out = process_file(f, out_dir, fps=args.fps, verbose=args.verbose, skip_existing=args.skip_existing)
                if out is not None:
                    print('Saved MP4:', out)
                else:
                    print('Skipped (already exists):', f)
            except Exception as e:
                print('Failed:', f, e)
    else:
        out = args.output if args.output is not None else None
        out_dir = os.path.dirname(out) if out is not None and os.path.dirname(out) != '' else None
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        try:
            outp = process_file(args.input, out, fps=args.fps, verbose=args.verbose, skip_existing=args.skip_existing)
            if outp is not None:
                print('Saved MP4:', outp)
            else:
                print('Skipped (already exists):', args.input)
        except Exception as e:
            print('Failed:', args.input, e)

if __name__ == '__main__':
    main()
