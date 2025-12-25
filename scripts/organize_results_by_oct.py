#!/usr/bin/env python3
"""
脚本名: organize_results_by_oct.py
功能: 将不同参数配置的OCT结果按tiff类型重新组织
作者: GitHub Copilot
日期: 2025年9月26日
说明:
    1. 遍历Output目录下的所有OCT数据文件夹
    2. 对每个OCT数据，在不同配置中查找对应的tiff文件
    3. 以OCT数据名创建文件夹，在其中按tiff类型创建子文件夹
    4. 每个子文件夹包含该OCT在不同配置下的同类型文件
"""

import os
import shutil
from pathlib import Path
from collections import defaultdict
import time

def fix_long_path(path):
    """修复Windows长路径问题"""
    path_str = str(path)
    if not path_str.startswith('\\\\?\\') and len(path_str) > 260:
        return Path('\\\\?\\' + path_str)
    return path

def get_short_config_name(config_name):
    """将长的配置名转换为短缩写"""
    # 移除 'ssdopu-' 前缀
    if config_name.startswith('ssdopu-'):
        short_name = config_name[7:]  # 移除 'ssdopu-'
    else:
        short_name = config_name

    # 创建缩写映射
    abbreviations = {
        'h1_5x5-sigma_3': 'h5x5s3',
        'h1_7x7-sigma_4': 'h7x7s4',
        'h1_9x9-sigma_5': 'h9x9s5',
        'kRL_03-kRU_21': 'k03021',
        'kRL_05-kRU_21': 'k05021',
        'kRL_07-kRU_21': 'k07021',
        'kRL_09-kRU_21': 'k09021'
    }

    # 如果是已知配置，使用缩写
    if short_name in abbreviations:
        return abbreviations[short_name]

    # 否则返回原名（但截断到合理长度）
    return short_name[:15] if len(short_name) > 15 else short_name

def main():
    # 配置参数
    output_root = Path(r'D:\1-Liu Jian\yongxin.wang\Output')
    organized_root = Path(r'D:\1-Liu Jian\yongxin.wang\Organized_Results')

    # 选择要组织的tiff文件类型 (可以根据需要修改)
    # 常见类型包括: 'Struc', 'dopu_SS', 'cumLA', 'drLA', 'PhR', 'Stokes' 等
    selected_types = [
        '1-1_Struc_flat',           # 结构图像(平坦化)
        '1-1_Struc',                # 结构图像
        '1-3_1rep-Stokes',          # 1次重复Stokes参数
        '1-3_4rep-Stokes',          # 4次重复Stokes参数
        '1-4_dopu_SS',              # DOPU图像
        '2-1_cumLA',                # 累积线性退偏
        '2-2_cumLA_hsvColoring',    # 累积线性退偏(HSV着色)
        '2-3_drLA',                 # 动态范围线性退偏
        '2-4_drLA_hsvColoring',     # 动态范围线性退偏(HSV着色)
        '2-5_PhR',                  # 相位retardation
        '2-6_cumLA_rmBG',           # 累积线性退偏(去背景)
        '2-7_cumLA_rmBG_hsvColoring', # 累积线性退偏(去背景, HSV着色)
        '2-8_drLA_rmBG',            # 动态范围线性退偏(去背景)
        '2-9_drLA_rmBG_hsvColoring', # 动态范围线性退偏(去背景, HSV着色)
        '2-10_PhR_rmBG'             # 相位retardation(去背景)
    ]

    # 预运行检查
    print('开始预运行检查...')

    # 检查输出根目录
    if not output_root.exists():
        raise FileNotFoundError(f'输出根目录不存在: {output_root}')

    # 创建组织结果根目录（确保父目录存在）
    try:
        organized_root.mkdir(parents=True, exist_ok=True)
        print(f'组织结果根目录就绪: {organized_root}')
    except Exception as e:
        print(f'警告: 无法创建组织结果根目录 {organized_root}: {e}')
        # 尝试使用当前目录的子目录
        organized_root = Path.cwd() / 'Organized_Results'
        organized_root.mkdir(parents=True, exist_ok=True)
        print(f'改用当前目录: {organized_root}')

    print('预运行检查完成！\n')

    # 获取所有参数配置文件夹
    print('扫描参数配置文件夹...')
    config_folders = [f for f in output_root.iterdir() if f.is_dir() and f.name.startswith('ssdopu-')]

    print(f'找到 {len(config_folders)} 个参数配置文件夹:')
    for i, folder in enumerate(config_folders, 1):
        print(f'  {i}. {folder.name}')
    print()

    # 收集所有OCT数据名称
    print('收集所有OCT数据名称...')
    oct_names = set()

    for config_folder in config_folders:
        oct_folders = [f for f in config_folder.iterdir() if f.is_dir()]
        for oct_folder in oct_folders:
            oct_names.add(oct_folder.name)

    oct_names = sorted(list(oct_names))
    print(f'找到 {len(oct_names)} 个OCT数据:')
    for i, oct_name in enumerate(oct_names, 1):
        print(f'  {i}. {oct_name}')
    print()

    # 组织结果
    print('开始组织结果...')
    total_start_time = time.time()

    organized_count = 0

    for oct_idx, oct_name in enumerate(oct_names, 1):
        print(f'组织OCT数据 {oct_idx}/{len(oct_names)}: {oct_name}')

        # 创建OCT结果文件夹
        oct_result_path = organized_root / oct_name
        try:
            oct_result_path.mkdir(parents=True, exist_ok=True)
            print(f'  创建OCT文件夹: {oct_result_path}')
        except Exception as e:
            print(f'  警告: 无法创建OCT文件夹 {oct_result_path}: {e}')
            continue

        # 为该OCT数据收集所有配置中的对应文件
        type_files_map = defaultdict(list)  # tiff类型 -> 文件列表

        for config_folder in config_folders:
            config_name = config_folder.name
            oct_folder_path = config_folder / oct_name
            tiff_path = oct_folder_path / 'tiff'

            # 检查tiff文件夹是否存在
            if not fix_long_path(tiff_path).exists():
                continue

            # 获取tiff文件
            tiff_files = list(fix_long_path(tiff_path).glob('*.tiff'))

            for tiff_file in tiff_files:
                file_name = tiff_file.name

                # 确定文件类型 - 使用改进的匹配逻辑
                file_type = None
                for selected_type in selected_types:
                    # 改进的匹配逻辑：处理不同的文件名模式
                    if selected_type.endswith('_hsvColoring'):
                        # 对于hsvColoring类型，查找 base_type-cfg*-hsvColoring_frame
                        base_type = selected_type.replace('_hsvColoring', '')
                        if f'_{base_type}-cfg' in file_name and 'hsvColoring_frame' in file_name:
                            file_type = selected_type
                            break
                    else:
                        # 对于普通类型，查找 _type-cfg 或 _type_frame
                        if f'_{selected_type}-cfg' in file_name or f'_{selected_type}_frame' in file_name:
                            file_type = selected_type
                            break

                # 如果没找到，尝试简单包含匹配（向后兼容）
                if file_type is None:
                    for selected_type in selected_types:
                        if selected_type in file_name:
                            file_type = selected_type
                            break

                if file_type:
                    # 存储该类型下的文件
                    type_files_map[file_type].append((tiff_file, config_name))

        # 为每个类型创建子文件夹并复制文件
        for tiff_type in selected_types:  # 改为遍历所有selected_types
            # 创建类型子文件夹
            type_result_path = oct_result_path / tiff_type
            try:
                type_result_path.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                print(f'    警告: 无法创建类型文件夹 {type_result_path}: {e}')
                continue

            # 获取该类型的所有文件
            file_list = type_files_map.get(tiff_type, [])

            # 复制文件
            copied_count = 0
            for tiff_file, config_name in file_list:
                file_name = tiff_file.name

                # 生成新文件名: 配置_类型信息
                short_config = get_short_config_name(config_name)
                
                # 从原文件名中提取类型信息（移除重复的OCT数据名称）
                # 文件名格式通常为: OCT名称_OCT名称_类型信息.tiff
                if file_name.count('_') >= 2:
                    # 去掉扩展名
                    name_without_ext = file_name.rsplit('.', 1)[0]
                    parts = name_without_ext.split('_')
                    
                    # 简化策略：找到第一个以1-或2-开头的部分，从那里开始提取类型信息
                    type_start_idx = -1
                    for i in range(len(parts)):
                        if parts[i].startswith(('1-', '2-')):
                            type_start_idx = i
                            break
                    
                    if type_start_idx >= 0:
                        # 提取从类型信息开始的所有部分
                        type_parts = parts[type_start_idx:]
                        type_info = '_'.join(type_parts)
                        new_filename = f"{short_config}_{type_info}.tiff"
                    else:
                        # 如果找不到1-或2-开头的部分，使用原文件名
                        new_filename = f"{short_config}_{file_name}"
                else:
                    # 简单文件名，直接使用
                    new_filename = f"{short_config}_{file_name}"
                
                dest_path = fix_long_path(type_result_path / new_filename)

                try:
                    # 检查源文件是否存在
                    if not fix_long_path(tiff_file).exists():
                        print(f'      警告: 源文件不存在 {tiff_file}')
                        continue

                    # 复制文件
                    shutil.copy2(str(fix_long_path(tiff_file)), str(fix_long_path(dest_path)))
                    copied_count += 1
                except Exception as e:
                    print(f'      警告: 复制文件失败 {tiff_file} -> {dest_path}: {e}')

            if copied_count > 0:
                print(f'    {tiff_type}: 复制了 {copied_count} 个文件')
            else:
                print(f'    {tiff_type}: 无文件')

        organized_count += 1

    # 完成总结
    total_time = time.time() - total_start_time
    print(f'\n{"="*40}')
    print('结果组织完成！')
    print(f'总耗时: {total_time:.1f}秒 ({total_time/60:.1f}分钟)')
    print(f'组织了 {organized_count} 个OCT数据')
    print(f'结果保存在: {organized_root}')
    print('='*40)

    # 生成使用指南
    print('\n使用指南:')
    print('1. 每个OCT数据都有独立的文件夹')
    print('2. 每个OCT文件夹下按tiff类型创建子文件夹')
    print('3. 每个类型子文件夹包含该OCT在不同配置下的同类型文件')
    print('4. 文件名格式: 配置名_原文件名')
    print('5. 支持的tiff类型: 结构图像、DOPU、Stokes参数、线性退偏、相位retardation等')
    print('6. 可以比较同一OCT数据在不同配置下的处理效果')

if __name__ == '__main__':
    main()