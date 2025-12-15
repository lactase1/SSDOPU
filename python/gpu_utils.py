"""
GPU加速工具模块
提供GPU/CPU自动切换功能，使用CuPy进行GPU加速
"""

import numpy as np
from typing import Union, Tuple, Optional
import warnings

# 尝试导入CuPy
try:
    import cupy as cp
    from cupyx.scipy import ndimage as cp_ndimage
    from cupyx.scipy.signal import fftconvolve as cp_fftconvolve
    GPU_AVAILABLE = True
    print("CuPy已加载，GPU加速可用")
except ImportError:
    GPU_AVAILABLE = False
    cp = None
    cp_ndimage = None
    cp_fftconvolve = None
    warnings.warn("CuPy未安装，将使用CPU计算。安装方法: pip install cupy-cuda12x")

# 选择计算后端
xp = cp if GPU_AVAILABLE else np


def get_array_module(arr):
    """
    获取数组对应的计算模块(numpy或cupy)
    
    Args:
        arr: 输入数组
    
    Returns:
        对应的计算模块
    """
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        return cp
    return np


def to_gpu(arr: np.ndarray) -> Union[np.ndarray, 'cp.ndarray']:
    """
    将numpy数组转移到GPU
    
    Args:
        arr: numpy数组
    
    Returns:
        GPU数组(如果GPU可用)或原数组
    """
    if GPU_AVAILABLE and isinstance(arr, np.ndarray):
        return cp.asarray(arr)
    return arr


def to_cpu(arr) -> np.ndarray:
    """
    将GPU数组转移到CPU
    
    Args:
        arr: GPU数组或numpy数组
    
    Returns:
        numpy数组
    """
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)
    return arr


def fft(arr, n: Optional[int] = None, axis: int = -1):
    """
    GPU加速的FFT
    
    Args:
        arr: 输入数组
        n: FFT长度
        axis: 计算轴
    
    Returns:
        FFT结果
    """
    xp_local = get_array_module(arr)
    return xp_local.fft.fft(arr, n=n, axis=axis)


def ifft(arr, n: Optional[int] = None, axis: int = -1):
    """
    GPU加速的IFFT
    
    Args:
        arr: 输入数组
        n: IFFT长度
        axis: 计算轴
    
    Returns:
        IFFT结果
    """
    xp_local = get_array_module(arr)
    return xp_local.fft.ifft(arr, n=n, axis=axis)


def hilbert(arr, axis: int = 0):
    """
    GPU加速的Hilbert变换
    
    Args:
        arr: 输入数组
        axis: 计算轴
    
    Returns:
        解析信号
    """
    xp_local = get_array_module(arr)
    N = arr.shape[axis]
    
    # FFT
    Xf = xp_local.fft.fft(arr, axis=axis)
    
    # 构建Hilbert滤波器
    h = xp_local.zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2
    
    # 扩展维度以便广播
    shape = [1] * arr.ndim
    shape[axis] = N
    h = h.reshape(shape)
    
    # 应用滤波器并IFFT
    return xp_local.fft.ifft(Xf * h, axis=axis)


def imfilter(arr, kernel, mode: str = 'reflect'):
    """
    GPU加速的图像滤波
    
    Args:
        arr: 输入图像
        kernel: 滤波核
        mode: 边界模式
    
    Returns:
        滤波结果
    """
    xp_local = get_array_module(arr)
    
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        kernel_gpu = cp.asarray(kernel)
        return cp_ndimage.convolve(arr, kernel_gpu, mode=mode)
    else:
        from scipy import ndimage
        return ndimage.convolve(arr, kernel, mode=mode)


def medfilt1(arr, kernel_size: int):
    """
    GPU加速的一维中值滤波
    
    Args:
        arr: 输入一维数组
        kernel_size: 滤波核大小
    
    Returns:
        滤波结果
    """
    xp_local = get_array_module(arr)
    
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        return cp_ndimage.median_filter(arr, size=kernel_size)
    else:
        from scipy.ndimage import median_filter
        return median_filter(arr, size=kernel_size)


def gaussian_filter(arr, sigma: float):
    """
    GPU加速的高斯滤波
    
    Args:
        arr: 输入数组
        sigma: 高斯标准差
    
    Returns:
        滤波结果
    """
    xp_local = get_array_module(arr)
    
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        return cp_ndimage.gaussian_filter(arr, sigma=sigma)
    else:
        from scipy.ndimage import gaussian_filter as scipy_gaussian
        return scipy_gaussian(arr, sigma=sigma)


def binary_dilation(arr, structure=None, iterations: int = 1):
    """
    GPU加速的二值膨胀
    
    Args:
        arr: 二值图像
        structure: 结构元素
        iterations: 迭代次数
    
    Returns:
        膨胀结果
    """
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        if structure is not None:
            structure = cp.asarray(structure)
        return cp_ndimage.binary_dilation(arr, structure=structure, iterations=iterations)
    else:
        from scipy.ndimage import binary_dilation as scipy_binary_dilation
        return scipy_binary_dilation(arr, structure=structure, iterations=iterations)


def binary_erosion(arr, structure=None, iterations: int = 1):
    """
    GPU加速的二值腐蚀
    
    Args:
        arr: 二值图像
        structure: 结构元素
        iterations: 迭代次数
    
    Returns:
        腐蚀结果
    """
    if GPU_AVAILABLE and isinstance(arr, cp.ndarray):
        if structure is not None:
            structure = cp.asarray(structure)
        return cp_ndimage.binary_erosion(arr, structure=structure, iterations=iterations)
    else:
        from scipy.ndimage import binary_erosion as scipy_binary_erosion
        return scipy_binary_erosion(arr, structure=structure, iterations=iterations)


def create_disk_structuring_element(radius: int):
    """
    创建圆形结构元素
    
    Args:
        radius: 半径
    
    Returns:
        圆形结构元素
    """
    size = 2 * radius + 1
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    return mask.astype(np.uint8)


def tukey_window(M: int, alpha: float = 0.5):
    """
    生成Tukey窗(余弦窗)
    
    Args:
        M: 窗口长度
        alpha: 形状参数 (0=矩形, 1=Hann窗)
    
    Returns:
        窗口数组
    """
    if alpha <= 0:
        return np.ones(M)
    elif alpha >= 1:
        return np.hanning(M)
    
    n = np.arange(M)
    width = int(np.floor(alpha * (M - 1) / 2.0))
    
    # 左边缘
    n1 = n[0:width + 1]
    window = np.zeros(M)
    window[0:width + 1] = 0.5 * (1 + np.cos(np.pi * (-1 + 2.0 * n1 / alpha / (M - 1))))
    
    # 中间平坦部分
    window[width + 1:M - width - 1] = 1
    
    # 右边缘
    n2 = n[M - width - 1:]
    window[M - width - 1:] = 0.5 * (1 + np.cos(np.pi * (-2.0 / alpha + 1 + 2.0 * n2 / alpha / (M - 1))))
    
    return window


def smooth1d(arr, window_size: int):
    """
    一维平滑滤波
    
    Args:
        arr: 输入一维数组
        window_size: 窗口大小
    
    Returns:
        平滑结果
    """
    xp_local = get_array_module(arr)
    kernel = xp_local.ones(window_size) / window_size
    return xp_local.convolve(arr, kernel, mode='same')


class GPUContext:
    """
    GPU上下文管理器，用于自动管理GPU内存
    """
    def __init__(self, device_id: int = 0):
        self.device_id = device_id
        
    def __enter__(self):
        if GPU_AVAILABLE:
            cp.cuda.Device(self.device_id).use()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if GPU_AVAILABLE:
            # 清理GPU内存
            cp.get_default_memory_pool().free_all_blocks()
        return False


def get_gpu_memory_info() -> dict:
    """
    获取GPU内存信息
    
    Returns:
        包含GPU内存信息的字典
    """
    if not GPU_AVAILABLE:
        return {"available": False}
    
    mempool = cp.get_default_memory_pool()
    return {
        "available": True,
        "used_bytes": mempool.used_bytes(),
        "total_bytes": mempool.total_bytes(),
        "n_free_blocks": mempool.n_free_blocks()
    }


if __name__ == "__main__":
    # 测试GPU功能
    print(f"GPU可用: {GPU_AVAILABLE}")
    
    if GPU_AVAILABLE:
        print(f"GPU内存信息: {get_gpu_memory_info()}")
        
        # 测试GPU计算
        with GPUContext():
            arr_cpu = np.random.rand(1000, 1000).astype(np.float32)
            arr_gpu = to_gpu(arr_cpu)
            
            # 测试FFT
            result_gpu = fft(arr_gpu, axis=0)
            result_cpu = to_cpu(result_gpu)
            
            print(f"FFT测试完成，结果形状: {result_cpu.shape}")
