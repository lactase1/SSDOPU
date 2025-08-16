import os
import numpy as np
import pydicom
from PIL import Image
import glob
import sys

def extract_frame_from_dcm(input_folder, frame_number=301):
    """
    从指定文件夹及其子文件夹中的所有DCM文件中提取指定帧并保存为TIFF文件
    
    参数:
    input_folder: 包含DCM文件的根文件夹路径
    frame_number: 要提取的帧号（从1开始计数），默认为301
    """
    
    # 检查输入文件夹是否存在
    if not os.path.exists(input_folder):
        print(f"错误: 指定的文件夹不存在: {input_folder}")
        return
    
    # 递归查找所有包含DCM文件的文件夹
    dcm_folders = []
    
    # 遍历所有子文件夹
    for root, dirs, files in os.walk(input_folder):
        # 检查当前文件夹是否包含DCM文件
        dcm_files = glob.glob(os.path.join(root, "*.dcm"))
        if dcm_files:
            dcm_folders.append((root, dcm_files))
    
    if not dcm_folders:
        print("未找到包含DCM文件的文件夹")
        return
    
    print(f"找到 {len(dcm_folders)} 个包含DCM文件的文件夹")
    
    # 处理每个包含DCM文件的文件夹
    for folder_path, dcm_files in dcm_folders:
        print(f"\n正在处理文件夹: {folder_path}")
        print(f"该文件夹包含 {len(dcm_files)} 个DCM文件")
        
        # 处理该文件夹中的每个DCM文件
        for dcm_file in dcm_files:
            try:
                # 读取DCM文件
                ds = pydicom.dcmread(dcm_file)
                
                # 获取像素数据
                pixel_data = ds.pixel_array
                
                # 检查帧数
                if len(pixel_data.shape) < 3:
                    print(f"警告: {os.path.basename(dcm_file)} 不包含多帧数据，跳过")
                    continue
                    
                total_frames = pixel_data.shape[0]
                print(f"文件 {os.path.basename(dcm_file)} 包含 {total_frames} 帧")
                
                # 检查帧号是否有效
                if frame_number > total_frames:
                    print(f"警告: {os.path.basename(dcm_file)} 只有 {total_frames} 帧，无法提取第 {frame_number} 帧，跳过")
                    continue
                
                # 提取指定帧（注意索引从0开始）
                frame_data = pixel_data[frame_number - 1]
                
                # 生成输出文件名
                base_name = os.path.splitext(os.path.basename(dcm_file))[0]
                output_filename = f"{base_name}_frame{frame_number}.tiff"
                output_path = os.path.join(folder_path, output_filename)  # 保存在同一个文件夹中
                
                # 保存为TIFF文件
                if frame_data.dtype == np.uint8:
                    image = Image.fromarray(frame_data)
                else:
                    # 对于其他数据类型，可能需要归一化
                    # 将数据归一化到0-255范围
                    frame_normalized = ((frame_data - frame_data.min()) / 
                                      (frame_data.max() - frame_data.min()) * 255).astype(np.uint8)
                    image = Image.fromarray(frame_normalized)
                
                image.save(output_path)
                print(f"已保存: {output_filename}")
                
            except Exception as e:
                print(f"处理 {os.path.basename(dcm_file)} 时出错: {str(e)}")
    
    print("\n所有文件处理完成!")

def main():
    if len(sys.argv) > 1:
        input_folder = sys.argv[1]
    else:
        input_folder = input("请输入包含DCM文件的根文件夹路径: ").strip()
    
    # 检查文件夹是否存在
    if not os.path.exists(input_folder):
        print(f"错误: 指定的文件夹不存在: {input_folder}")
        return
    
    # 提取第301帧
    extract_frame_from_dcm(input_folder, frame_number=301)

# 使用示例
if __name__ == "__main__":
    main()