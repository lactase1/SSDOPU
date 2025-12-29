# MATLAB Coder 输入类型整理

本文档整理了 `function/` 目录下 MATLAB 函数的输入类型，用于 MATLAB Coder 生成 C 代码。每个函数包含函数签名、输入描述和 `coder.typeof` 模板。

## 通用说明
- 将模板中的占位符（如 `*_max`）替换为实际的最大维度值（整数）。
- 所有类型假设使用 `double`（实数或复数），除非特别注明。
- 对于结构体，使用 `coder.typeof(struct(...), [1 1], [0 0])` 定义固定字段。
- 运行 `codegen` 时，使用 `-config:lib` 生成静态库。

## 1. calculateSplitSpectrumDOPU
**函数签名：**
```matlab
[dopu_splitSpectrum, dopu_ss] = calculateSplitSpectrumDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg)
```

**输入类型：**
- `Bd1`, `Bd2`：复数 double，尺寸 [SPL × nX × nr]
  ```matlab
  Bd_type = coder.typeof(1+1i, [SPL_max, nX_max, nr_max], [1 1 1]);
  ```
- `params`：结构体，包含 `dopu.do_ssdopu`、`dopu.verbose`、`processing.show_img` 等布尔字段
  ```matlab
  params_type = coder.typeof(struct('dopu', struct('do_ssdopu', false, 'verbose', false), 'processing', struct('show_img', false)), [1 1], [0 0]);
  ```
- `SPL, nX, nr, nWin, winL`：标量 int32
  ```matlab
  SPL_type = coder.typeof(int32(0), [1 1], [0 0]);
  nX_type  = coder.typeof(int32(0), [1 1], [0 0]);
  nr_type  = coder.typeof(int32(0), [1 1], [0 0]);
  nWin_type = coder.typeof(int32(0), [1 1], [0 0]);
  winL_type = coder.typeof(int32(0), [1 1], [0 0]);
  ```
- `windex`：int32 向量，尺寸 [1 × nWin]
  ```matlab
  windex_type = coder.typeof(int32(0), [1, nWin_max], [0, 1]);
  ```
- `winG`：double 向量，尺寸 [winL × 1]
  ```matlab
  winG_type = coder.typeof(1.0, [winL_max, 1], [0 0]);
  ```
- `czrg`：int32 向量，尺寸 [nZcrop × 1]（nZcrop ≤ SPL）
  ```matlab
  czrg_type = coder.typeof(int32(0), [SPL_max, 1], [1 0]);
  ```

**codegen 示例：**
```matlab
codegen calculateSplitSpectrumDOPU -args {Bd_type, Bd_type, params_type, SPL_type, nX_type, nr_type, nWin_type, windex_type, winL_type, winG_type, czrg_type} -config:lib
```

## 2. calculateCombinedDOPU
**函数签名：**
```matlab
dopu_combined = calculateCombinedDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg)
```

**输入类型：**
- 同 `calculateSplitSpectrumDOPU`，但 `params` 需包含 `dopu.do_combined`。

## 3. calculateSpatialDOPU
**函数签名：**
```matlab
dopu_spatial = calculateSpatialDOPU(IMG_ch1, IMG_ch2, params)
```

**输入类型：**
- `IMG_ch1`, `IMG_ch2`：复数 double，尺寸 [nZ × nX × nRep]
  ```matlab
  IMG_type = coder.typeof(1+1i, [nZ_max, nX_max, nRep_max], [1 1 1]);
  ```
- `params`：结构体，包含 `dopu.do_spatial` 等布尔字段
  ```matlab
  params_type = coder.typeof(struct('dopu', struct('do_spatial', false)), [1 1], [0 0]);
  ```

**codegen 示例：**
```matlab
codegen calculateSpatialDOPU -args {IMG_type, IMG_type, params_type} -config:lib
```

## 4. cumulativeQUV
**函数签名：**
```matlab
[S0, S1, S2, S3] = cumulativeQUV(IMG_ch1, IMG_ch2)
```

**输入类型：**
- `IMG_ch1`, `IMG_ch2`：复数 double，尺寸 [nZ × nX] 或 [nZ × nX × nRep]
  ```matlab
  IMG_type = coder.typeof(1+1i, [nZ_max, nX_max, nRep_max], [1 1 1]);
  ```

**codegen 示例：**
```matlab
codegen cumulativeQUV -args {IMG_type, IMG_type} -config:lib
```

## 5. FreeSpace_PSOCT_Optimized
**函数签名：**
```matlab
[LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw] = FreeSpace_PSOCT_Optimized(Qm, Um, Vm, Stru, test_seg_top, h1, h2, Avnum, dopuMap, enableDopuPhaseSupp)
```

**输入类型：**
- `Qm, Um, Vm, Stru`：double，尺寸 [nZ × nX]
  ```matlab
  M_type = coder.typeof(1.0, [nZ_max, nX_max], [1 1]);
  ```
- `test_seg_top`：int32 向量，尺寸 [1 × nX]
  ```matlab
  test_seg_type = coder.typeof(int32(0), [1, nX_max], [0 1]);
  ```
- `h1, h2`：double 矩阵，尺寸 [h_size × h_size]（小矩阵，如 5×5）
  ```matlab
  h_type = coder.typeof(1.0, [h_max, h_max], [0 0]);
  ```
- `Avnum`：标量 int32
  ```matlab
  Avnum_type = coder.typeof(int32(0), [1 1], [0 0]);
  ```
- `dopuMap`：double 矩阵，尺寸 [nZ × nX]（可为空）
  ```matlab
  dopuMap_type = coder.typeof(1.0, [nZ_max, nX_max], [1 1]);
  ```
- `enableDopuPhaseSupp`：logical 标量
  ```matlab
  en_type = coder.typeof(true, [1 1], [0 0]);
  ```

**codegen 示例：**
```matlab
codegen FreeSpace_PSOCT_Optimized -args {M_type, M_type, M_type, M_type, test_seg_type, h_type, h_type, Avnum_type, dopuMap_type, en_type} -config:lib
```

## 6. vWinAvgFiltOpt
**函数签名：**
```matlab
outFrameWfilt = vWinAvgFiltOpt(inFrame, inWeight, kRL, kRU, Nsec)
```

**输入类型：**
- `inFrame, inWeight`：double，尺寸 [nZ × nX]
  ```matlab
  frame_type = coder.typeof(1.0, [nZ_max, nX_max], [1 1]);
  ```
- `kRL, kRU, Nsec`：标量 int32
  ```matlab
  k_type = coder.typeof(int32(0), [1 1], [0 0]);
  ```

**codegen 示例：**
```matlab
codegen vWinAvgFiltOpt -args {frame_type, frame_type, k_type, k_type, k_type} -config:lib
```

**注意：** 该函数使用 `fspecial` 和动态 cell，可能需修改为固定滤波核数组以兼容 MATLAB Coder。

## 7. overwrite_config_params
**函数签名：**
```matlab
overwrite_config_params(updates, config_file)
```

**输入类型：**
- `updates`：结构体，包含更新字段
  ```matlab
  updates_type = coder.typeof(struct('polarization', struct('Avnum', int32(0)), 'filters', struct('h2_sigma', 1.0, 'h2_size', [int32(0), int32(0)])), [1 1], [0 0]);
  ```
- `config_file`：字符串
  ```matlab
  file_type = coder.typeof('string', [1, 256], [0 1]);
  ```

**注意：** 该函数涉及文件 I/O，不适合直接生成纯 C 库。建议跳过或重写为纯字符串处理。

## 示例 Wrapper 脚本
创建一个 MATLAB 脚本（如 `codegen_wrapper.m`），定义所有类型变量，然后运行 `codegen`。

```matlab
% 示例：为 calculateSplitSpectrumDOPU 生成 C 代码
SPL_max = 1024; nX_max = 512; nr_max = 3; nWin_max = 9; winL_max = 256;
Bd_type = coder.typeof(1+1i, [SPL_max, nX_max, nr_max], [1 1 1]);
params_type = coder.typeof(struct('dopu', struct('do_ssdopu', false, 'verbose', false), 'processing', struct('show_img', false)), [1 1], [0 0]);
SPL_type = coder.typeof(int32(0), [1 1], [0 0]);
nX_type = coder.typeof(int32(0), [1 1], [0 0]);
nr_type = coder.typeof(int32(0), [1 1], [0 0]);
nWin_type = coder.typeof(int32(0), [1 1], [0 0]);
windex_type = coder.typeof(int32(0), [1, nWin_max], [0, 1]);
winL_type = coder.typeof(int32(0), [1 1], [0 0]);
winG_type = coder.typeof(1.0, [winL_max, 1], [0 0]);
czrg_type = coder.typeof(int32(0), [SPL_max, 1], [1 0]);

codegen calculateSplitSpectrumDOPU -args {Bd_type, Bd_type, params_type, SPL_type, nX_type, nr_type, nWin_type, windex_type, winL_type, winG_type, czrg_type} -config:lib
```

将 `_max` 替换为实际值。