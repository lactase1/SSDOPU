"""
PS-OCT处理主入口脚本
批量处理OCT文件
"""

import os
import sys
import time
from pathlib import Path
from typing import List, Optional

from config_params import config_params
from psoct_processor import process_single_file


def main(data_path: Optional[str] = None, output_base: Optional[str] = None):
    """
    主处理函数
    
    Args:
        data_path: 数据目录路径
        output_base: 输出目录路径
    """
    # 默认路径
    if data_path is None:
        data_path = r'D:\1-Liu Jian\yongxin.wang\PSOCT'
    if output_base is None:
        output_base = r'D:\1-Liu Jian\yongxin.wang\tmp'
    
    data_path = Path(data_path)
    output_base = Path(output_base)
    
    # 检查数据路径
    if not data_path.exists():
        raise FileNotFoundError(f"数据路径不存在: {data_path}")
    
    # 创建输出路径
    output_base.mkdir(parents=True, exist_ok=True)
    print(f"输出目录: {output_base}")
    
    # 获取所有OCT文件
    print("正在处理数据...")
    oct_files = list(data_path.glob("*.oct"))
    
    if not oct_files:
        print(f"在 {data_path} 中未找到OCT文件")
        return
    
    # 显示找到的文件
    print(f"找到 {len(oct_files)} 个 OCT 文件:")
    for i, f in enumerate(oct_files, 1):
        print(f"[{i}] {f.name}")
    
    # 加载配置参数
    params = config_params()
    params.print_summary()
    
    # 逐个处理文件
    total_start_time = time.time()
    success_count = 0
    
    for i, oct_file in enumerate(oct_files, 1):
        print(f"\n{'='*50}")
        print(f"开始处理第 {i}/{len(oct_files)} 个文件: {oct_file.name}")
        print(f"{'='*50}")
        
        try:
            start_time = time.time()
            
            # 调用单文件处理函数
            process_single_file(str(oct_file), str(output_base), params)
            
            proc_time = time.time() - start_time
            print(f"文件 {oct_file.name} 处理成功, 耗时: {proc_time:.2f} 秒 ({proc_time/60:.2f} 分钟)")
            success_count += 1
            
        except Exception as e:
            print(f"处理文件 {oct_file.name} 时出错: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # 总结
    total_time = time.time() - total_start_time
    print(f"\n{'='*50}")
    print(f"所有文件处理完成!")
    print(f"成功: {success_count}/{len(oct_files)}")
    print(f"总耗时: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")
    print(f"{'='*50}")


if __name__ == "__main__":
    # 从命令行参数获取路径
    if len(sys.argv) >= 3:
        data_path = sys.argv[1]
        output_base = sys.argv[2]
    elif len(sys.argv) == 2:
        data_path = sys.argv[1]
        output_base = None
    else:
        data_path = None
        output_base = None
    
    main(data_path, output_base)
