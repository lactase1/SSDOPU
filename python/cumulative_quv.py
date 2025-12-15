"""
Stokes参数计算模块
从两通道复数OCT图像计算Stokes参数(S0, S1, S2, S3)
"""

import numpy as np
from typing import Tuple, Union

try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None

from gpu_utils import get_array_module, to_gpu, to_cpu


def cumulative_quv(img_ch1: np.ndarray, img_ch2: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    从两通道复数OCT图像计算Stokes参数
    
    算法原理:
        - Stokes参数是描述偏振光的完整参数集
        - S0: 总光强度 |E1|² + |E2|²
        - S1: 水平/垂直偏振差 |E1|² - |E2|²
        - S2: ±45°偏振差 2|E1||E2|cos(θ)
        - S3: 左/右圆偏振差 2|E1||E2|sin(-θ)
        - θ为两通道信号的相位差
    
    Args:
        img_ch1: 通道1的复数OCT信号矩阵 [Z×X或Z×X×Rep]
        img_ch2: 通道2的复数OCT信号矩阵 [Z×X或Z×X×Rep]
    
    Returns:
        S0, S1, S2, S3: 四个Stokes参数
    """
    xp = get_array_module(img_ch1)
    
    # 计算两通道信号的相位差
    phase = xp.angle(img_ch2 * xp.conj(img_ch1))
    
    # 计算幅度
    abs_ch1 = xp.abs(img_ch1)
    abs_ch2 = xp.abs(img_ch2)
    
    # 计算幅度平方
    abs_ch1_sq = abs_ch1 ** 2
    abs_ch2_sq = abs_ch2 ** 2
    
    # 计算四个Stokes参数
    S0 = abs_ch1_sq + abs_ch2_sq                           # 总光强度
    S1 = abs_ch1_sq - abs_ch2_sq                           # 水平-垂直偏振分量差
    S2 = 2.0 * abs_ch1 * abs_ch2 * xp.cos(phase)           # +45°与-45°偏振分量差
    S3 = 2.0 * abs_ch1 * abs_ch2 * xp.sin(-phase)          # 右旋与左旋圆偏振分量差
    
    return S0, S1, S2, S3


def cumulative_quv_gpu(img_ch1: np.ndarray, img_ch2: np.ndarray) -> Tuple:
    """
    GPU加速版本的Stokes参数计算（返回GPU数组）

    Args:
        img_ch1: 通道1的复数OCT信号矩阵 (numpy或cupy数组)
        img_ch2: 通道2的复数OCT信号矩阵 (numpy或cupy数组)

    Returns:
        S0, S1, S2, S3: 四个Stokes参数（cupy数组）
    """
    if not GPU_AVAILABLE:
        return cumulative_quv(img_ch1, img_ch2)

    # 转移到GPU（如果尚未在GPU上）
    img_ch1_gpu = to_gpu(img_ch1)
    img_ch2_gpu = to_gpu(img_ch2)

    # 计算（返回GPU数组）
    S0, S1, S2, S3 = cumulative_quv(img_ch1_gpu, img_ch2_gpu)

    return S0, S1, S2, S3


if __name__ == "__main__":
    # 测试
    print("测试Stokes参数计算...")
    
    # 创建测试数据
    nZ, nX = 320, 500
    img_ch1 = np.random.randn(nZ, nX) + 1j * np.random.randn(nZ, nX)
    img_ch2 = np.random.randn(nZ, nX) + 1j * np.random.randn(nZ, nX)
    
    # CPU计算
    S0, S1, S2, S3 = cumulative_quv(img_ch1, img_ch2)
    print(f"CPU计算完成: S0形状={S0.shape}")
    
    # GPU计算
    if GPU_AVAILABLE:
        S0_gpu, S1_gpu, S2_gpu, S3_gpu = cumulative_quv_gpu(img_ch1, img_ch2)
        print(f"GPU计算完成: S0形状={S0_gpu.shape}")
        
        # 验证结果（将GPU结果转回CPU进行对比）
        diff = np.abs(S0 - cp.asnumpy(S0_gpu)).max()
        print(f"CPU/GPU差异: {diff}")
