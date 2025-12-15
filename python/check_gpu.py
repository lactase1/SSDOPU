"""
GPU加速检查和性能测试脚本
检查CuPy是否正确安装以及GPU加速效果
"""

import numpy as np
import time

print("="*60)
print("PS-OCT GPU加速状态检查")
print("="*60)

# 检查CuPy
print("\n[1] 检查CuPy安装...")
try:
    import cupy as cp
    print(f"✓ CuPy版本: {cp.__version__}")
    
    # 检查CUDA
    print(f"✓ CUDA版本: {cp.cuda.runtime.runtimeGetVersion()}")
    
    # 检查GPU设备
    device_count = cp.cuda.runtime.getDeviceCount()
    print(f"✓ 检测到 {device_count} 个GPU设备")
    
    for i in range(device_count):
        cp.cuda.Device(i).use()
        props = cp.cuda.runtime.getDeviceProperties(i)
        print(f"  - GPU {i}: {props['name'].decode()}")
        
        # 显示显存
        mem_info = cp.cuda.Device(i).mem_info
        free_mem = mem_info[0] / 1024**3
        total_mem = mem_info[1] / 1024**3
        print(f"    显存: {free_mem:.2f}GB 可用 / {total_mem:.2f}GB 总计")
    
    GPU_AVAILABLE = True
except Exception as e:
    print(f"✗ CuPy未正确安装或GPU不可用: {e}")
    print("\n安装CuPy:")
    print("  - CUDA 11.x: pip install cupy-cuda11x")
    print("  - CUDA 12.x: pip install cupy-cuda12x")
    GPU_AVAILABLE = False

# 简单性能测试
if GPU_AVAILABLE:
    print("\n[2] GPU vs CPU 性能测试...")
    
    # 创建测试数据
    size = (1000, 1000)
    print(f"测试数据大小: {size}")
    
    # CPU测试: FFT
    print("\n测试1: FFT运算")
    data_cpu = np.random.randn(*size) + 1j * np.random.randn(*size)
    
    start = time.time()
    for _ in range(10):
        result_cpu = np.fft.fft2(data_cpu)
    cpu_time = time.time() - start
    print(f"  CPU: {cpu_time:.3f}秒 (10次)")
    
    # GPU测试: FFT
    data_gpu = cp.asarray(data_cpu)
    cp.cuda.Stream.null.synchronize()  # 确保数据传输完成
    
    start = time.time()
    for _ in range(10):
        result_gpu = cp.fft.fft2(data_gpu)
    cp.cuda.Stream.null.synchronize()  # 确保计算完成
    gpu_time = time.time() - start
    print(f"  GPU: {gpu_time:.3f}秒 (10次)")
    print(f"  加速比: {cpu_time/gpu_time:.1f}x")
    
    # 测试2: 矩阵运算
    print("\n测试2: 矩阵乘法")
    A_cpu = np.random.randn(*size)
    B_cpu = np.random.randn(*size)
    
    start = time.time()
    for _ in range(10):
        C_cpu = A_cpu @ B_cpu
    cpu_time = time.time() - start
    print(f"  CPU: {cpu_time:.3f}秒 (10次)")
    
    A_gpu = cp.asarray(A_cpu)
    B_gpu = cp.asarray(B_cpu)
    cp.cuda.Stream.null.synchronize()
    
    start = time.time()
    for _ in range(10):
        C_gpu = A_gpu @ B_gpu
    cp.cuda.Stream.null.synchronize()
    gpu_time = time.time() - start
    print(f"  GPU: {gpu_time:.3f}秒 (10次)")
    print(f"  加速比: {cpu_time/gpu_time:.1f}x")
    
    # 测试3: 高斯滤波
    print("\n测试3: 高斯滤波")
    from scipy.ndimage import gaussian_filter as scipy_gaussian
    from cupyx.scipy.ndimage import gaussian_filter as cupy_gaussian
    
    image_cpu = np.random.randn(2000, 2000)
    
    start = time.time()
    for _ in range(5):
        filtered_cpu = scipy_gaussian(image_cpu, sigma=5)
    cpu_time = time.time() - start
    print(f"  CPU: {cpu_time:.3f}秒 (5次)")
    
    image_gpu = cp.asarray(image_cpu)
    cp.cuda.Stream.null.synchronize()
    
    start = time.time()
    for _ in range(5):
        filtered_gpu = cupy_gaussian(image_gpu, sigma=5)
    cp.cuda.Stream.null.synchronize()
    gpu_time = time.time() - start
    print(f"  GPU: {gpu_time:.3f}秒 (5次)")
    print(f"  加速比: {cpu_time/gpu_time:.1f}x")

# 检查项目文件
print("\n[3] 检查项目GPU化状态...")
import os

files_to_check = [
    ("freespace_psoct_vectorized.py", "向量化DDG算法"),
    ("psoct_processor.py", "主处理模块"),
    ("v_win_avg_filt_opt.py", "自适应滤波"),
    ("gpu_utils.py", "GPU工具函数"),
]

for filename, desc in files_to_check:
    if os.path.exists(f"python/{filename}"):
        print(f"✓ {desc} ({filename})")
    else:
        print(f"✗ {desc} ({filename}) - 文件不存在")

print("\n[4] 使用建议...")
if GPU_AVAILABLE:
    print("✓ GPU加速已就绪！")
    print("\n运行处理:")
    print("  python main.py <数据目录> <输出目录>")
    print("\n监控GPU使用:")
    print("  nvidia-smi -l 1  # 实时监控GPU状态")
else:
    print("✗ GPU不可用，将使用CPU计算（速度较慢）")
    print("\n请安装CuPy以启用GPU加速:")
    print("  pip install cupy-cuda12x  # 或 cupy-cuda11x")

print("\n" + "="*60)
