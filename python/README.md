# PS-OCT Processing - Python Version (GPU Accelerated)

## 简介

这是PS-OCT（偏振敏感光学相干断层扫描）处理的Python版本，从MATLAB代码转换而来，支持GPU加速。

## 功能特点

- **GPU加速**: 使用CuPy实现CUDA GPU加速，大幅提升处理速度
- **OCT文件读取**: 支持自定义.oct文件格式
- **DOPU计算**: 
  - 分裂谱DOPU (Split-Spectrum DOPU)
  - 空间DOPU (Spatial DOPU)
  - 组合DOPU (Combined DOPU)
- **Stokes参数计算**: 计算完整的偏振Stokes参数 (S0, S1, S2, S3)
- **DDG算法**: 深度依赖梯度法计算光轴和相位延迟
- **表面分割**: 自动检测组织表面
- **DICOM输出**: 保存结果为DICOM格式

## 安装

### 1. 创建虚拟环境（推荐）

```bash
python -m venv venv
# Windows
.\venv\Scripts\activate
# Linux/Mac
source venv/bin/activate
```

### 2. 安装依赖

```bash
pip install -r requirements.txt
```

### 3. 安装GPU支持（可选但推荐）

根据你的CUDA版本选择合适的CuPy包：

```bash
# CUDA 11.x
pip install cupy-cuda11x

# CUDA 12.x
pip install cupy-cuda12x
```

## 使用方法

### 命令行

```bash
# 处理单个目录下的所有OCT文件
python main.py <数据目录> <输出目录>

# 示例
python main.py "D:\PSOCT\data" "D:\PSOCT\output"
```

### Python代码

```python
from main import main
from config_params import config_params

# 使用默认参数
main(data_path="D:/PSOCT/data", output_base="D:/PSOCT/output")

# 自定义参数
params = config_params()
params.processing.max_frames = 100  # 限制处理帧数
params.dopu.do_combined = True      # 使用组合DOPU

from psoct_processor import process_single_file
process_single_file("path/to/file.oct", "output/dir", params)
```

## 项目结构

```
python/
├── __init__.py              # 包初始化
├── main.py                  # 主入口脚本
├── config_params.py         # 配置参数
├── psoct_processor.py       # 主处理模块
├── oct_file_reader.py       # OCT文件读取
├── gpu_utils.py             # GPU工具函数
├── cumulative_quv.py        # Stokes参数计算
├── calculate_split_spectrum_dopu.py  # 分裂谱DOPU
├── calculate_spatial_dopu.py         # 空间DOPU
├── calculate_combined_dopu.py        # 组合DOPU
├── v_win_avg_filt_opt.py    # 自适应DOPU滤波
├── freespace_psoct.py       # DDG算法核心
├── image_utils.py           # 图像处理工具
├── dicom_utils.py           # DICOM处理
├── requirements.txt         # 依赖列表
└── README.md                # 本文件
```

## 配置参数

主要参数在 `config_params.py` 中定义：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `tiff.make_tiff` | True | 是否生成TIFF文件 |
| `tiff.save_dicom` | True | 是否保存DICOM文件 |
| `processing.max_frames` | 0 | 最大处理帧数 (0=全部) |
| `dopu.do_combined` | True | 使用组合DOPU |
| `dopu.do_ssdopu` | True | 启用分裂谱DOPU |
| `polarization.avnum` | 3 | DDG平均层数 |

## GPU加速

当安装了CuPy时，程序会自动使用GPU加速。可以通过以下方式检查GPU状态：

```python
from gpu_utils import GPU_AVAILABLE, get_gpu_memory_info

print(f"GPU可用: {GPU_AVAILABLE}")
if GPU_AVAILABLE:
    print(get_gpu_memory_info())
```

## 输出文件

处理完成后会在输出目录生成以下文件：

- `*_1-1_Struc.dcm` - 结构图像
- `*_1-4_dopu_SS.dcm` - DOPU图像
- `*_2-1_cumLA-*.dcm` - 累积双折射
- `*_2-5_PhR-*.dcm` - 相位延迟
- `tiff/` - TIFF格式的提取帧

## 与MATLAB版本的对应关系

| MATLAB函数 | Python模块 |
|------------|-----------|
| `config_params.m` | `config_params.py` |
| `cumulativeQUV.m` | `cumulative_quv.py` |
| `calculateSplitSpectrumDOPU.m` | `calculate_split_spectrum_dopu.py` |
| `calculateSpatialDOPU.m` | `calculate_spatial_dopu.py` |
| `calculateCombinedDOPU.m` | `calculate_combined_dopu.py` |
| `vWinAvgFiltOpt.m` | `v_win_avg_filt_opt.py` |
| `FreeSpace_PSOCT_3_DDG_rmBG_7.m` | `freespace_psoct.py` |

## 性能优化建议

1. **使用GPU**: 安装CuPy可显著提升处理速度
2. **限制帧数**: 测试时可设置 `max_frames` 参数
3. **调整worker数**: 可在 `params.parallel.max_workers` 中设置

## 许可证

内部使用

## 联系方式

如有问题，请联系开发团队。
