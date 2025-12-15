# PS-OCT Processing Python Package
# GPU加速的偏振敏感光学相干断层扫描处理

"""
PS-OCT处理Python包

主要功能:
- OCT文件读取 (.oct格式)
- DOPU (偏振均匀度) 计算
  - 分裂谱DOPU
  - 空间DOPU
  - 组合DOPU
- Stokes参数计算
- DDG算法 (深度依赖梯度)
- 表面分割
- 光学衰减系数计算
- DICOM文件读写
- GPU加速 (可选，使用CuPy)

使用方法:
    from python import main
    main.main(data_path, output_base)

或命令行:
    python main.py <数据目录> <输出目录>

依赖:
    必需: numpy, scipy, Pillow
    推荐: cupy (GPU加速), pydicom (DICOM支持)

作者: PS-OCT处理系统
版本: 1.0
"""

__version__ = "1.0.0"
__author__ = "PS-OCT Processing System"

from .config_params import config_params, ConfigParams
from .psoct_processor import process_single_file
from .main import main

__all__ = [
    'config_params',
    'ConfigParams',
    'process_single_file',
    'main',
]
