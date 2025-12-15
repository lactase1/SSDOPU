# PS-OCT GPU加速优化说明

## 🚀 优化总结

已完成对PS-OCT处理系统的GPU向量化优化，主要改进：

### 1. 核心算法完全向量化 ✅
- **文件**: `freespace_psoct_vectorized.py`
- **改进**: 移除了所有逐像素的Python循环
- **方法**: 
  - 使用CuPy数组切片实现滑动窗口
  - 批量计算旋转矩阵（Rodrigues公式）
  - 向量化DDG切向法相位计算
  - 并行化三点法主轴拟合
- **预期提速**: **10-100倍**

### 2. 数据流全程GPU化 ✅
- **文件**: `psoct_processor.py`
- **改进**: 
  - 数据读取后立即上传到GPU
  - FFT、Hilbert变换、Stokes计算全部在GPU上
  - 仅在最后保存结果时转回CPU
- **避免**: CPU-GPU间频繁数据传输

### 3. 辅助函数GPU化 ✅
- **自适应DOPU滤波** (`v_win_avg_filt_opt.py`): 
  - 接受cupy数组输入
  - 使用`cupyx.scipy.ndimage.gaussian_filter`
  - 输出保持在GPU
- **GPU工具函数** (`gpu_utils.py`):
  - 正确使用CuPy FFT
  - Hilbert变换基于CuPy实现
  - 自动检测数组类型切换计算后端

## 📦 使用方法

### 1. 检查GPU状态
```bash
cd python
python check_gpu.py
```

这会显示:
- CuPy和CUDA版本
- 可用GPU设备
- 显存信息
- GPU vs CPU性能对比

### 2. 运行处理
```bash
python main.py <OCT数据目录> <输出目录>
```

### 3. 监控GPU使用
在另一个终端运行:
```bash
nvidia-smi -l 1  # 每秒刷新
```

应该看到:
- **GPU利用率**: 80-100%（不再是0%）
- **显存使用**: 2-8GB（根据数据大小）
- **处理速度**: 比之前快10-100倍

## 🔧 优化细节

### 向量化策略

#### 原代码（慢）:
```python
for j in range(nX):  # 500次
    for i in range(output_depth):  # 300次
        # 单点计算
        P1 = planeData[0, :]
        P2 = planeData[Avnum-1, :]
        T1 = P2 - P1
        # ... 更多计算
```

#### 优化后（快）:
```python
# 一次性处理所有像素
P1 = ps_bg_rm[:, :, 0, :]  # [output_depth, nX, 3]
P2 = ps_bg_rm[:, :, Avnum-1, :]
T1 = P2 - P1  # 向量化减法
# GPU并行计算所有点
```

### 数据流优化

#### 原代码:
```
读取数据 → NumPy (CPU)
    ↓
处理 → NumPy (CPU)
    ↓
DDG算法 → 循环中频繁 to_gpu/to_cpu
    ↓
保存 → NumPy (CPU)
```

#### 优化后:
```
读取数据 → NumPy (CPU)
    ↓
立即上传 → CuPy (GPU) ← 只传输一次
    ↓
FFT → CuPy (GPU)
    ↓
Stokes → CuPy (GPU)
    ↓
DOPU → CuPy (GPU)
    ↓
DDG (向量化) → CuPy (GPU)
    ↓
最后下载 → NumPy (CPU) ← 只传输一次
    ↓
保存
```

## 📊 预期性能提升

| 操作 | 原速度 | 优化后速度 | 提速比 |
|------|-------|----------|--------|
| DDG算法 | 很慢（循环） | 快 | **50-100x** |
| FFT | 中等 | 快 | **5-10x** |
| 高斯滤波 | 中等 | 快 | **10-20x** |
| 整体处理 | 慢 | 快 | **预计10-30x** |

## ⚠️ 注意事项

### 显存管理
- 单个B-scan约需1-2GB显存
- 如果显存不足，会出现`Out of Memory`错误
- 解决方法:
  ```python
  # 在config_params.py中限制帧数
  params.processing.max_frames = 50  # 减少处理帧数
  ```

### CuPy版本
确保CuPy与CUDA版本匹配:
```bash
# 查看CUDA版本
nvidia-smi

# 安装对应的CuPy
pip install cupy-cuda12x  # CUDA 12.x
# 或
pip install cupy-cuda11x  # CUDA 11.x
```

### 调试
如果遇到错误:
1. 先运行 `python check_gpu.py` 确认GPU可用
2. 检查CUDA和CuPy版本是否匹配
3. 用小数据测试 (`max_frames=10`)
4. 查看nvidia-smi确认显存未满

## 📝 修改的文件列表

1. **新建**: `freespace_psoct_vectorized.py` - 完全向量化的DDG算法
2. **修改**: `psoct_processor.py` - GPU数据流优化
3. **修改**: `v_win_avg_filt_opt.py` - GPU版自适应滤波
4. **修改**: `gpu_utils.py` - 确认GPU函数正确性
5. **新建**: `check_gpu.py` - GPU状态检查工具
6. **新建**: `GPU_OPTIMIZATION_GUIDE.md` - 本文档

## 🎯 测试建议

1. **小规模测试**:
   ```python
   # 在main.py中临时添加
   params.processing.max_frames = 10
   ```

2. **对比测试**:
   - 运行优化前的代码记录时间
   - 运行优化后的代码记录时间
   - 对比GPU利用率

3. **验证结果**:
   - 确认输出图像与原MATLAB版本一致
   - 检查DICOM文件正确生成

## 💡 进一步优化建议

如果需要更高性能:
1. **Batch Processing**: 同时处理多个B-scan
2. **Multi-GPU**: 使用多张GPU卡并行
3. **混合精度**: 使用float16减少显存占用
4. **CUDA Kernel**: 为关键操作编写自定义CUDA kernel

## 📞 问题排查

### GPU利用率仍然很低?
- 检查是否使用了新的`freespace_psoct_vectorized`
- 确认CuPy正确安装: `python -c "import cupy; print(cupy.__version__)"`
- 查看是否有错误回退到CPU

### 显存不足?
- 减少`max_frames`
- 减小处理的Z范围
- 关闭其他占用GPU的程序

### 结果不一致?
- 检查数值精度（GPU可能有微小差异）
- 验证索引转换（0-based vs 1-based）
- 对比中间结果

---

**开发团队**: PS-OCT Processing System  
**优化日期**: 2025年12月15日  
**GPU加速框架**: CuPy
