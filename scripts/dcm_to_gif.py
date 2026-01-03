#!/usr/bin/env python3
"""
dcm_to_gif.py

Convert DICOM (.dcm) multi-frame files to animated GIF(s).

Usage:
    python dcm_to_gif.py <input_path> [--output <out_path_or_dir>] [--duration MS] [--loop N] [--verbose]

If <input_path> is a file, produces one GIF. If it's a directory, processes all .dcm files inside and
writes per-file GIFs to the output directory (or same folder if not given).

Dependencies: pydicom, Pillow, numpy

Examples:
    python dcm_to_gif.py data/example.dcm --output out.gif --duration 100
    python dcm_to_gif.py data/dcm_folder --output out_dir --duration 80 --loop 0
"""

import os
import argparse
import math
import numpy as np
from PIL import Image
import pydicom


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
    """Normalize numpy frame to uint8 2D or 3D array suitable for Image.fromarray"""
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


def frames_to_gif(frames, out_path, duration=100, loop=0, verbose=False):
    if len(frames) == 0:
        raise ValueError('No frames to save')
    pil_frames = []
    for i, f in enumerate(frames):
        f_u8 = normalize_frame_to_uint8(f)
        # If grayscale 2D -> convert to RGB for better GIF results
        if f_u8.ndim == 2:
            img = Image.fromarray(f_u8).convert('RGBA')
        elif f_u8.ndim == 3 and f_u8.shape[2] == 3:
            img = Image.fromarray(f_u8).convert('RGBA')
        elif f_u8.ndim == 3 and f_u8.shape[2] == 4:
            img = Image.fromarray(f_u8).convert('RGBA')
        else:
            # unexpected channels: reduce to first channel
            img = Image.fromarray(f_u8[..., 0]).convert('RGBA')
        # Convert to P (palette) to produce smaller GIFs
        pal = img.convert('P', palette=Image.ADAPTIVE)
        pil_frames.append(pal)
        if verbose:
            print(f'Prepared frame {i+1}/{len(frames)} size={img.size}')

    first = pil_frames[0]
    others = pil_frames[1:]
    # Save GIF
    first.save(out_path, format='GIF', save_all=True, append_images=others, duration=duration, loop=loop)
    if verbose:
        print(f'GIF saved: {out_path} ({len(frames)} frames, duration={duration}ms, loop={loop})')


def process_file(dcm_path, out_path=None, duration=100, loop=0, verbose=False):
    ds = pydicom.dcmread(dcm_path)
    frames = get_frames_from_dataset(ds, verbose=verbose)
    if out_path is None:
        out_path = os.path.splitext(dcm_path)[0] + '.gif'
    elif os.path.isdir(out_path):
        out_path = os.path.join(out_path, os.path.splitext(os.path.basename(dcm_path))[0] + '.gif')

    frames_to_gif(frames, out_path, duration=duration, loop=loop, verbose=verbose)
    return out_path


def main():
    parser = argparse.ArgumentParser(description='Convert DICOM (multi-frame) to animated GIF')
    parser.add_argument('input', help='DICOM file or directory containing .dcm files')
    parser.add_argument('--output', '-o', help='Output GIF file path or directory (defaults to input folder)', default=None)
    parser.add_argument('--duration', '-d', type=int, default=100, help='Frame duration in milliseconds')
    parser.add_argument('--loop', '-l', type=int, default=0, help='Number of GIF loops (0=infinite)')
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
                out = process_file(f, out_dir, duration=args.duration, loop=args.loop, verbose=args.verbose)
                print('Saved GIF:', out)
            except Exception as e:
                print('Failed:', f, e)
    else:
        out = args.output if args.output is not None else None
        out_dir = os.path.dirname(out) if out is not None and os.path.dirname(out) != '' else None
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        try:
            outp = process_file(args.input, out, duration=args.duration, loop=args.loop, verbose=args.verbose)
            print('Saved GIF:', outp)
        except Exception as e:
            print('Failed:', args.input, e)

if __name__ == '__main__':
    main()
