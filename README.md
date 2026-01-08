# SSDOPU_PSOCT_Code

![MATLAB](https://img.shields.io/badge/MATLAB-R2020b+-blue.svg)
![Python](https://img.shields.io/badge/Python-3.7+-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

这是一个基于 MATLAB 的PSOCT数据处理平台，能处理OCT数据，生成结构图像，偏振特性图像，可输出dcm格式，tiff格式，png格式的结果。


## 📑 目录

- [核心功能](#-核心功能)
- [项目结构](#-项目结构)
- [环境要求](#-环境要求)
- [快速开始](#-快速开始)
- [配置说明](#-配置说明)
- [高级功能](#-高级功能)
- [辅助工具](#-辅助工具)
- [输出说明](#-输出说明)

---

## 🌟 核心功能

### 偏振分析
- **多种 DOPU 计算**：
  - **分光谱 DOPU** (Split-spectrum)：支持9段半重叠或5段无重叠模式
  - **空间 DOPU** (Spatial)：基于3×3邻域平均计算
  - **组合 DOPU** (Combined)：结合分光谱和空间方法
- **双折射分析**：累积和局部光轴计算，相位延迟计算
- **偏振特性颜色图**：偏振特性彩色编码输出

### 图像处理
- **高性能处理**：支持分批流式处理（Batch Streaming），大幅降低内存占用
- **自适应滤波**：基于分裂谱DOPU质量的智能滤波
- **智能展平系统**：
  - 自动生成边界：没有边界会自动计算生成
  - 加载边界信息：基于 `topLines` 对齐数据
  - 后巩膜展平：支持外部边界MAT文件精准对齐
- **并行计算**：支持多核并行处理，显著提升处理速度

### 数据输出
- **多格式支持**：通过配置能自动生成 DCM（3D体数据及 En-face 切片）和 TIFF 图像
- **批处理模式**：支持一次性处理多个 OCT 文件
- **可视化**：自动生成结构图、DOPU图、双折射图等多种可视化结果

---

## 📁 项目结构

```
SSDOPU/
├── README.md                          # 项目说明文档
├── split_dopu_psoct_main.m           # 主处理脚本
│
├── function/                          # 核心处理函数库
│   ├── config_params.m               # 配置参数文件（重要）
│   ├── calculateSplitSpectrumDOPU.m  # 分光谱DOPU计算
│   ├── calculateSpatialDOPU.m        # 空间DOPU计算
│   ├── calculateCombinedDOPU.m       # 组合DOPU计算
│   ├── calLAPhRALL.m                 # 双折射和相位延迟计算
│   ├── calOAC.m                      # 光轴计算
│   ├── cumulativeQUV.m               # 累积偏振态计算
│   ├── filterVectorFieldAdaptive.m   # 自适应向量场滤波
│   ├── vWinAvgFiltOpt.m              # 优化的窗口平均滤波
│   ├── surf_seg.m                    # 表面分割
│   ├── FreeSpace_PSOCT_Optimized.m   # PS-OCT优化处理
│   ├── convert_dcm_to_tiff_local.m   # DCM转TIFF转换
│   ├── overwrite_config_params.m     # 配置参数覆写
│   ├── print_progress.m              # 进度条显示
│   └── quColoring.m                  # 偏振态着色
│
└── scripts/                           # 辅助脚本工具
    ├── organize_results_by_oct.py    # 结果文件整理工具
    ├── dcm_to_gif.py                 # DCM转GIF动画
    ├── data_helper.py                # 数据处理辅助函数
    ├── tiff.py                       # TIFF处理工具
    ├── extract_dcm_frame.m           # DCM帧提取
    ├── generate_structure_pngs.m     # 批量生成结构图PNG
    ├── parameter_sweep_script.m      # 参数扫描脚本
    └── run_planefit_sweep.m          # 平面拟合参数扫描
```

### 文件夹功能说明

#### `function/` - 核心函数库
包含所有PS-OCT数据处理的核心算法和函数，包括：
- DOPU计算算法（分光谱、空间、组合）
- 偏振分析函数（双折射、相位延迟、光轴）
- 图像处理函数（滤波、分割、展平）
- 配置管理和工具函数

#### `scripts/` - 辅助工具集
提供数据管理、格式转换、批处理等辅助功能：
- Python脚本：用于结果整理、格式转换、数据可视化
- MATLAB脚本：用于参数扫描、批量处理、结果生成

---

## 💻 环境要求

### MATLAB环境

### Python环境（用于辅助脚本）

## 🚀 快速开始

### 1. 环境配置
下载好matlab就行

### 2. 设置输入输出路径

编辑 `split_dopu_psoct_main.m` 文件：

```matlab
% ========== 设置数据路径 ==========
% 输入路径：包含.oct文件的文件夹
data_path = 'C:\your\data\path\';

% 输出路径：处理结果保存位置
output_base = 'D:\your\output\path\';
```

### 3. 配置处理参数

编辑 `function/config_params.m` 文件，调整关键参数（详见[配置说明](#-配置说明)）。

### 4. 运行处理

在MATLAB命令行中运行：
```matlab
run('split_dopu_psoct_main.m')
```

或直接打开 `split_dopu_psoct_main.m` 并点击运行按钮。

---

## ⚙️ 配置说明

所有处理参数集中在 `function/config_params.m` 文件中。以下是关键配置项说明：

### 基础处理参数

```matlab
%% ========== 基础处理参数 ==========
params.processing.disp_coef = -20.1;           % 色散补偿系数
params.processing.do_PhComp = 1;               % 是否进行相位/色散补偿 (1:启用, 0:禁用)
params.processing.do_medianshift = 1;          % 是否进行中值偏移校正
params.processing.do_reg = 0;                  % 是否进行帧间配准(运动去除)
params.processing.max_frames = 0;              % 最大处理帧数 (0:处理所有帧)
```

**说明**：
- `disp_coef`：色散补偿系数，根据OCT系统特性调整
- `do_PhComp`：启用相位补偿可提高图像质量，但会增加处理时间
- `max_frames`：用于快速测试，设为0处理全部帧

### DOPU计算配置

```matlab
%% ========== 分光谱DOPU设置 ==========
params.dopu.do_ssdopu = 1;                     % 是否启用分光谱DOPU (1:启用, 0:禁用)
params.dopu.ss_mode = 'overlap9';              % 分裂谱模式: 'overlap9' 或 'nonoverlap5'

%% ========== 空间DOPU设置 ==========
params.dopu.do_spatial = 0;                    % 是否启用空间DOPU

%% ========== 组合DOPU设置 ==========
params.dopu.do_combined = 1;                   % 是否启用组合DOPU
```

**说明**：
- `do_ssdopu`：分光谱DOPU是核心功能，建议启用
- `ss_mode`：`'overlap9'` 提供更好的噪声抑制，`'nonoverlap5'` 速度更快
- `do_combined`：组合方法通常提供最佳结果

### 偏振分析参数

```matlab
%% ========== 偏振分析参数 ==========
params.polarization.Avnum = 3;                 % 平均层数 (影响空间分辨率和信噪比)
params.polarization.kRL_cfg1 = 3;              % 滤波核下限
params.polarization.kRU_cfg1 = 21;             % 滤波核上限
params.polarization.PRRrg = [0.06 0.35];       % 相位延迟范围
```

**说明**：
- `Avnum`：增大可提高信噪比，但降低轴向分辨率（典型值：1-5）
- `kRL_cfg1` / `kRU_cfg1`：自适应滤波的核大小范围，影响滤波强度
- `PRRrg`：相位延迟显示范围，根据样品特性调整

### 滤波器设置

```matlab
%% ========== 高斯滤波核设置 ==========
params.filters.h1_size = [5 5];                % 小核尺寸 (细节保持)
params.filters.h1_sigma = 1.5;                 % 小核标准差

params.filters.h2_size = [19 19];              % 大核尺寸 (背景平滑)
params.filters.h2_sigma = 3;                   % 大核标准差
```

**说明**：
- 小核（h1）：用于保持细节特征
- 大核（h2）：用于抑制深层噪声
- 巩膜成像建议：h1保持较小（5×5），h2适当增大（15×15至21×21）

### 并行处理设置

```matlab
%% ========== 并行处理设置 ==========
params.parallel.maxWorkers = 47;               % 最大并行worker数
params.parallel.maxMemGB = 100;                % 最大可用内存(GB)
params.parallel.autoClosePool = false;         % 是否自动关闭并行池
```

**说明**：
- `maxWorkers`：设为0自动检测，或手动指定核心数
- `maxMemGB`：根据系统内存容量设置
- `autoClosePool`：设为false可在多次运行间保持并行池，节省启动时间

### 输出控制

```matlab
%% ========== TIFF生成控制参数 ==========
params.tiff.make_tiff = 0;                     % 是否生成TIFF文件
params.tiff.tiff_frame = 35;                   % 要提取的帧号
params.tiff.saveDicom = 1;                     % 是否保存DICOM文件

params.processing.enable_flatten_enface = 0;   % 是否生成展平En-face图像
params.processing.enable_enface_noflat = 0;    % 是否生成非展平En-face图像
```

---

## 🔧 高级功能

### 后巩膜边界展平功能

该功能允许基于外部提供的边界数据对齐巩膜结构，用于精确的巩膜层分析。

#### 1. 准备边界文件
准备一个 `.mat` 文件，包含与 B-scan 对应的边界位置数据（单位：像素）。
- **默认变量名**：`bottomLines`
- **维度**：应与 A-lines 数量及 B-scan 数量匹配（例如 `500×500` 或 `500×30`）

#### 2. 配置参数
在 `function/config_params.m` 中设置：
```matlab
%% ========== 后巩膜边界处理参数 ==========
params.files.sclera_boundary_path = 'C:\path\to\sclera_boundary.mat';
params.files.sclera_boundary_var = '';  % 留空自动探测，或指定变量名
```

#### 3. 输出说明
程序将自动生成以 `_sclera_flat` 结尾的文件：
- `*_4-1_Enface_Struc_sclera_flat.dcm`：基于后巩膜对齐的结构 En-face
- `*_4-2_Enface_cumLA_sclera_flat.dcm`：基于后巩膜对齐的双折射 En-face

### 自适应滤波策略

```matlab
%% ========== 输出端自适应滤波参数 ==========
params.filters.enable_output_adaptive = 0;     % 启用输出端混合滤波
params.filters.output_dopu_threshold = 0.35;   % DOPU阈值
params.filters.kRL_output = 13;                % 自适应滤波核下限
params.filters.kRU_output = 31;                % 自适应滤波核上限
```

**工作原理**：
- DOPU ≥ 阈值（高质量区域）：使用固定h2滤波，保留细节
- DOPU < 阈值（低质量区域）：使用自适应滤波，增强降噪

---

## 🛠️ 辅助工具

### Python脚本

#### `organize_results_by_oct.py`
将不同参数配置的处理结果按OCT文件重新组织。

**使用方法**：
```bash
python scripts/organize_results_by_oct.py
```

**功能**：
- 遍历输出目录中的所有配置文件夹
- 按OCT文件名创建文件夹
- 按TIFF类型分类组织结果

#### `dcm_to_gif.py`
将多帧DICOM文件转换为动画GIF。

**使用方法**：
```bash
# 转换单个文件
python scripts/dcm_to_gif.py input.dcm --output output.gif --duration 100

# 批量转换文件夹
python scripts/dcm_to_gif.py input_folder --output output_folder --duration 80
```

**参数说明**：
- `--duration`：帧间延时（毫秒）
- `--loop`：循环次数（0表示无限循环）
- `--verbose`：显示详细信息

### MATLAB脚本

#### `parameter_sweep_script.m`
参数扫描脚本，用于批量测试不同参数组合。

**功能**：
- 自动修改 `config_params.m` 中的参数
- 批量运行主处理脚本
- 保存不同参数配置的结果到独立文件夹

#### `generate_structure_pngs.m`
批量生成结构图PNG文件。

**功能**：
- 从DCM文件提取特定帧
- 生成高质量PNG图像
- 支持批处理多个文件

---

## 📊 输出说明

### 输出目录结构

```
output_base/
├── run_20260108_143022/              # 批处理时间戳目录
│   ├── sample001_processed/          # 单个OCT文件处理结果
│   │   ├── 1-Tiff_frame35/          # TIFF文件
│   │   ├── 2-Volume_3D/             # 3D体数据DCM
│   │   ├── 3-Enface/                # En-face切片
│   │   ├── 4-Flatten_Enface/        # 展平En-face
│   │   └── processing_log.txt       # 处理日志
│   └── sample002_processed/
└── ...
```

### 主要输出文件

#### TIFF文件 (`1-Tiff_frame35/`)
- `*_Structure.tif`：结构图
- `*_DOPU.tif`：DOPU图
- `*_cumLA.tif`：累积双折射图
- `*_PhR.tif`：相位延迟图

#### 3D体数据 (`2-Volume_3D/`)
- `*_1-1_Volume_Struc.dcm`：结构3D体
- `*_1-2_Volume_cumLA.dcm`：双折射3D体
- `*_1-3_Volume_DOPU.dcm`：DOPU 3D体

#### En-face切片 (`3-Enface/`)
- `*_2-1_Enface_Struc.dcm`：结构En-face
- `*_2-2_Enface_cumLA.dcm`：双折射En-face
- `*_2-3_Enface_DOPU.dcm`：DOPU En-face

#### 展平数据 (`4-Flatten_Enface/`)
- `*_4-1_Enface_Struc_flat.dcm`：展平结构En-face
- `*_4-2_Enface_cumLA_flat.dcm`：展平双折射En-face
- `*_4-1_Enface_Struc_sclera_flat.dcm`：后巩膜展平（如果启用）

---

## 📝 常见问题

### Q: 内存不足怎么办？
A: 调整以下参数：
- 减小 `params.parallel.maxWorkers`
- 增大 `params.parallel.batchSize`（分批处理）
- 设置 `params.processing.max_frames` 限制处理帧数

### Q: 处理速度太慢？
A: 尝试以下优化：
- 启用并行处理（设置 `maxWorkers > 0`）
- 禁用不需要的输出（`make_tiff = 0`）
- 选择更快的DOPU模式（`ss_mode = 'nonoverlap5'`）

### Q: 图像质量不理想？
A: 调整以下参数：
- 增大 `Avnum` 提高信噪比
- 调整滤波核大小（`h1_size`, `h2_size`）
- 修改 `kRL_cfg1` 和 `kRU_cfg1` 范围

---

## 📄 许可证

MIT License

---

## 👥 贡献者

- 王永鑫 (Yongxin Wang)

---

## 📮 联系方式

如有问题或建议，请通过以下方式联系：
- 提交 Issue
- 发送邮件

---

**最后更新时间**：2026年1月8日
- `*_4-3_Enface_PhR_sclera_flat.dcm`: 基于后巩膜对齐的延迟相位 En-face。

---

## 📁 项目结构

- `split_dopu_psoct_main.m`: 主程序入口，负责批处理流逻辑。
- `function/`: 核心算法库（DOPU、滤波、彩色编码、分割等）。
- `python/`: 相关的 Python 处理工具和 GPU 优化模块。
- `scripts/`: 数据整理与参数搜索辅助脚本。

## ⚠️ 注意事项

- **内存优化**：本项目默认使用 `single` 精度处理数据，若处理超大规模数据，请适当调小 `params.parallel.batchSize`。
- **并行计算**：程序会自动启动并行池（Parpool），工作线程数由 `params.parallel.maxWorkers` 控制。
