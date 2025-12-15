"""
DICOM文件处理模块
读写DICOM格式的医学图像
"""

import numpy as np
from pathlib import Path
from typing import Union, Optional, List
import warnings

try:
    import pydicom
    from pydicom.dataset import Dataset, FileDataset
    from pydicom.uid import generate_uid
    PYDICOM_AVAILABLE = True
except ImportError:
    PYDICOM_AVAILABLE = False
    warnings.warn("pydicom未安装，DICOM功能不可用。安装方法: pip install pydicom")


def dicom_write(
    data: np.ndarray,
    filename: str,
    bits_allocated: int = 8,
    photometric: str = 'MONOCHROME2'
):
    """
    将numpy数组保存为DICOM文件
    
    Args:
        data: 图像数据，支持以下格式:
            - 2D: [H, W] 单帧灰度图
            - 3D: [H, W, C] 单帧RGB图 或 [H, W, N] 多帧灰度图
            - 4D: [H, W, C, N] 多帧RGB图
        filename: 输出文件名
        bits_allocated: 位深度 (8 或 16)
        photometric: 光度解释 ('MONOCHROME2' 或 'RGB')
    """
    if not PYDICOM_AVAILABLE:
        raise ImportError("pydicom未安装")
    
    # 确保数据类型正确
    if bits_allocated == 8:
        data = data.astype(np.uint8)
    elif bits_allocated == 16:
        data = data.astype(np.uint16)
    
    # 创建文件元信息
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.7'  # Secondary Capture
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    
    # 创建数据集
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble=b"\0" * 128)
    
    # 设置基本属性
    ds.SOPClassUID = file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.Modality = 'OT'  # Other
    ds.PatientName = 'Anonymous'
    ds.PatientID = '000000'
    
    # 根据数据维度设置属性
    if data.ndim == 2:
        # 单帧灰度图
        ds.Rows, ds.Columns = data.shape
        ds.NumberOfFrames = 1
        ds.PhotometricInterpretation = photometric
        ds.SamplesPerPixel = 1
        ds.BitsAllocated = bits_allocated
        ds.BitsStored = bits_allocated
        ds.HighBit = bits_allocated - 1
        ds.PixelRepresentation = 0
        ds.PixelData = data.tobytes()
        
    elif data.ndim == 3:
        if data.shape[2] == 3:
            # 单帧RGB图
            ds.Rows, ds.Columns = data.shape[:2]
            ds.NumberOfFrames = 1
            ds.PhotometricInterpretation = 'RGB'
            ds.SamplesPerPixel = 3
            ds.PlanarConfiguration = 0
            ds.BitsAllocated = bits_allocated
            ds.BitsStored = bits_allocated
            ds.HighBit = bits_allocated - 1
            ds.PixelRepresentation = 0
            ds.PixelData = data.tobytes()
        else:
            # 多帧灰度图
            ds.Rows, ds.Columns = data.shape[:2]
            ds.NumberOfFrames = data.shape[2]
            ds.PhotometricInterpretation = photometric
            ds.SamplesPerPixel = 1
            ds.BitsAllocated = bits_allocated
            ds.BitsStored = bits_allocated
            ds.HighBit = bits_allocated - 1
            ds.PixelRepresentation = 0
            # 重排数据: [H, W, N] -> [N, H, W]
            data_reordered = np.transpose(data, (2, 0, 1))
            ds.PixelData = data_reordered.tobytes()
            
    elif data.ndim == 4:
        # 多帧RGB图
        ds.Rows, ds.Columns = data.shape[:2]
        ds.NumberOfFrames = data.shape[3]
        ds.PhotometricInterpretation = 'RGB'
        ds.SamplesPerPixel = 3
        ds.PlanarConfiguration = 0
        ds.BitsAllocated = bits_allocated
        ds.BitsStored = bits_allocated
        ds.HighBit = bits_allocated - 1
        ds.PixelRepresentation = 0
        # 重排数据: [H, W, C, N] -> [N, H, W, C]
        data_reordered = np.transpose(data, (3, 0, 1, 2))
        ds.PixelData = data_reordered.tobytes()
    
    # 保存文件
    ds.save_as(filename)


def dicom_read(filename: str) -> np.ndarray:
    """
    读取DICOM文件为numpy数组
    
    Args:
        filename: DICOM文件路径
    
    Returns:
        图像数据数组
    """
    if not PYDICOM_AVAILABLE:
        raise ImportError("pydicom未安装")
    
    ds = pydicom.dcmread(filename)
    return ds.pixel_array


def convert_dcm_to_tiff(
    dcm_folder: str,
    target_frame: int,
    output_prefix: str = ''
):
    """
    将DCM文件转换为TIFF文件
    
    Args:
        dcm_folder: DCM文件所在文件夹
        target_frame: 目标帧号 (1-based)
        output_prefix: 输出文件前缀
    """
    from PIL import Image
    
    if not PYDICOM_AVAILABLE:
        raise ImportError("pydicom未安装")
    
    dcm_folder = Path(dcm_folder)
    if not dcm_folder.exists():
        raise FileNotFoundError(f"DCM文件夹不存在: {dcm_folder}")
    
    # 创建输出目录
    tiff_output_dir = dcm_folder / 'tiff'
    tiff_output_dir.mkdir(exist_ok=True)
    
    # 查找所有DCM文件
    dcm_files = list(dcm_folder.glob('*.dcm'))
    
    if not dcm_files:
        warnings.warn(f"在指定目录中未找到DCM文件: {dcm_folder}")
        return
    
    print(f"转换 {len(dcm_files)} 个DCM文件的第{target_frame}帧到TIFF...")
    
    success_count = 0
    for dcm_file in dcm_files:
        try:
            # 读取DCM文件
            dcm_data = dicom_read(str(dcm_file))
            
            # 确定帧的位置
            if dcm_data.ndim == 4:
                total_frames = dcm_data.shape[0]
                frame_dim = 0
            elif dcm_data.ndim == 3:
                if dcm_data.shape[2] == 3:
                    # RGB图像
                    frame_dim = None
                    total_frames = 1
                else:
                    total_frames = dcm_data.shape[0]
                    frame_dim = 0
            else:
                frame_dim = None
                total_frames = 1
            
            # 提取目标帧
            if frame_dim is None:
                frame_data = dcm_data
            else:
                actual_frame = min(target_frame - 1, total_frames - 1)  # 转为0索引
                frame_data = dcm_data[actual_frame]
            
            # 生成输出文件名
            base_name = dcm_file.stem
            if output_prefix:
                tiff_filename = f"{output_prefix}_{base_name}_frame{target_frame}.tiff"
            else:
                tiff_filename = f"{base_name}_frame{target_frame}.tiff"
            
            tiff_filepath = tiff_output_dir / tiff_filename
            
            # 保存为TIFF
            img = Image.fromarray(frame_data)
            img.save(str(tiff_filepath))
            
            success_count += 1
            
        except Exception as e:
            print(f"处理文件 {dcm_file.name} 出错: {e}")
            continue
    
    print(f"完成! 成功转换 {success_count}/{len(dcm_files)} 个文件到: {tiff_output_dir}")


if __name__ == "__main__":
    print("测试DICOM模块...")
    
    if PYDICOM_AVAILABLE:
        # 测试写入
        test_data = np.random.randint(0, 255, (100, 100), dtype=np.uint8)
        dicom_write(test_data, 'test_output.dcm')
        print("DICOM写入测试完成")
        
        # 测试读取
        read_data = dicom_read('test_output.dcm')
        print(f"DICOM读取测试完成，数据形状: {read_data.shape}")
        
        # 清理测试文件
        Path('test_output.dcm').unlink()
    else:
        print("pydicom未安装，跳过测试")
