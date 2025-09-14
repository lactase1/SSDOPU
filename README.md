# PSOCT 数据处理 GPU 加速版本

## 概述
本项目包含一组 MATLAB 函数和脚本，用于处理 PSOCT (偏振敏感光学相干断层成像) 数据。这些函数已经优化为可以利用 NVIDIA GPU 进行加速计算，特别适合 RTX 4090 等高性能 GPU。

## 主要功能
- 使用 GPU 加速处理大量 PSOCT 数据
- 支持分段光谱 (split-spectrum) DOPU 计算
- 实现了两种配置的 LA 和 Phase 计算
- 使用可变高斯滤波器优化 Q U V 处理
- 支持自动表面分割和背景移除
- 生成和保存 DICOM 格式的结果

## 文件结构
- `rPSOCT_06_splitSpec_3_3D_gpu.m` - 主处理函数，GPU 加速版本
- `process_PSOCT_with_GPU.m` - 批处理脚本，用于处理多个文件
- `function/` 目录 - 包含辅助函数:
  - `cumulativeQUV_gpu.m` - GPU 加速的 Stokes 参数计算
  - `vWinAvgFiltOpt_2_1_gpu.m` - GPU 加速的可变窗口高斯滤波
  - `calLAPhRALL_gpu.m` - GPU 加速的线偏振和相位视网膜计算 (配置1)
  - `calLAPhRcfg2_gpu.m` - GPU 加速的线偏振和相位视网膜计算 (配置2)
  - 其他辅助函数...

## 使用方法

### 系统要求
- MATLAB R2020b 或更高版本
- Parallel Computing Toolbox
- CUDA 兼容的 NVIDIA GPU (推荐 RTX 4090)
- 至少 16GB GPU 内存 (处理大型数据集)

### 快速开始
1. 打开 MATLAB
2. 设置工作目录为包含脚本的文件夹
3. 运行 `process_PSOCT_with_GPU.m` 脚本
4. 根据提示选择数据文件夹和输出文件夹
5. 等待处理完成

### 参数调整
在 `process_PSOCT_with_GPU.m` 中，可以根据需要调整以下参数:
- `do_avg` - 是否执行平均处理
- `do_eig` - 是否执行特征值处理
- `do_cfg1` - 是否使用配置1
- `do_cfg2` - 是否使用配置2
- `do_PhComp` - 是否执行相位补偿
- `do_ssdopu` - 是否使用 split-spectrum DOPU
- `show_img` - 是否显示处理过程中的图像
- `saveDicom` - 是否保存 DICOM 格式的结果

### 处理单个文件
也可以直接调用主函数处理单个文件:
```matlab
rPSOCT_06_splitSpec_3_3D_gpu(...
    'file', '数据文件路径.oct', ...
    'outDir', '输出目录', ...
    'do_avg', true, ...
    'do_cfg1', true, ...
    'do_ssdopu', true);
```

## 性能建议
- 对于 RTX 4090，推荐使用 CUDA 12.0 或更高版本
- 确保系统有足够的 RAM (推荐 32GB+)
- 使用 SSD 存储以提高数据读写速度
- 关闭其他 GPU 密集型应用程序以获得最佳性能

## 故障排除
- 如果遇到 "GPU 内存不足" 错误，尝试减小处理的数据尺寸
- 如果未检测到 GPU，确保已安装正确的 CUDA 驱动程序和 Parallel Computing Toolbox
- 对于大文件，可能需要增加 MATLAB 的 JVM 堆内存

## 作者
[Yongxin Wang]