"""
空间DOPU计算模块
基于3x3邻域平均计算DOPU
"""

import numpy as np
from typing import Optional

try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None

from gpu_utils import get_array_module, to_gpu, to_cpu
from cumulative_quv import cumulative_quv


def calculate_spatial_dopu(
    img_ch1: np.ndarray,
    img_ch2: np.ndarray,
    params,
    use_gpu: bool = True
) -> np.ndarray:
    """
    计算空间DOPU (基于3x3邻域平均)
    
    算法原理:
        - 计算每个像素的Stokes参数 (S0, S1, S2, S3)
        - 对每个深度层，使用3x3邻域对Stokes参数进行空间平均
        - 从平均的Stokes参数计算DOPU
        - DOPU = sqrt((<S1>/<S0>)² + (<S2>/<S0>)² + (<S3>/<S0>)²)
    
    Args:
        img_ch1: 通道1的复数OCT信号矩阵 [Z×X×Rep]
        img_ch2: 通道2的复数OCT信号矩阵 [Z×X×Rep]
        params: 参数结构体
        use_gpu: 是否使用GPU加速
    
    Returns:
        dopu_spatial: 空间DOPU结果 [Z×X]
    """
    # 参数验证
    if img_ch1.shape != img_ch2.shape:
        raise ValueError("img_ch1和img_ch2的尺寸必须相同")
    
    # 检查是否需要进行空间DOPU计算
    if not params.dopu.do_spatial:
        nZ, nX = img_ch1.shape[:2]
        return np.ones((nZ, nX))
    
    # 获取数据维度
    if img_ch1.ndim == 3:
        nZ, nX, nRep = img_ch1.shape
    else:
        nZ, nX = img_ch1.shape
        nRep = 1
    
    # 选择计算后端
    if use_gpu and GPU_AVAILABLE:
        xp = cp
        img_ch1 = to_gpu(img_ch1)
        img_ch2 = to_gpu(img_ch2)
    else:
        xp = np
    
    # 对重复次数进行平均 (如果有多个重复)
    if nRep > 1:
        img_ch1 = xp.mean(img_ch1, axis=2)
        img_ch2 = xp.mean(img_ch2, axis=2)
    
    # 计算Stokes参数
    S0, S1, S2, S3 = cumulative_quv(img_ch1, img_ch2)
    
    # 初始化空间DOPU结果
    dopu_spatial = xp.zeros((nZ, nX))
    
    # 3点均值核
    kernel = xp.ones(3) / 3
    
    # 对每个深度层进行空间滤波
    for z in range(nZ):
        # 获取当前层的Stokes参数
        S0_layer = S0[z, :]
        S1_layer = S1[z, :]
        S2_layer = S2[z, :]
        S3_layer = S3[z, :]
        
        # 边界扩充处理
        S0_padded = xp.pad(S0_layer, 1, mode='edge')
        S1_padded = xp.pad(S1_layer, 1, mode='edge')
        S2_padded = xp.pad(S2_layer, 1, mode='edge')
        S3_padded = xp.pad(S3_layer, 1, mode='edge')
        
        # 使用3点均值滤波器进行空间平均
        S0_avg = xp.convolve(S0_padded, kernel, mode='valid')
        S1_avg = xp.convolve(S1_padded, kernel, mode='valid')
        S2_avg = xp.convolve(S2_padded, kernel, mode='valid')
        S3_avg = xp.convolve(S3_padded, kernel, mode='valid')
        
        # 避免除零错误
        S0_avg_safe = xp.maximum(S0_avg, xp.finfo(float).eps)
        
        # 计算DOPU
        dopu_spatial[z, :] = xp.sqrt(
            (S1_avg / S0_avg_safe)**2 +
            (S2_avg / S0_avg_safe)**2 +
            (S3_avg / S0_avg_safe)**2
        )
    
    # 转回CPU
    if use_gpu and GPU_AVAILABLE:
        dopu_spatial = to_cpu(dopu_spatial)
    
    return dopu_spatial


if __name__ == "__main__":
    from config_params import config_params
    
    print("测试空间DOPU计算...")
    
    # 创建测试参数
    params = config_params()
    params.dopu.do_spatial = True
    
    # 创建测试数据
    nZ, nX, nRep = 320, 500, 4
    img_ch1 = np.random.randn(nZ, nX, nRep) + 1j * np.random.randn(nZ, nX, nRep)
    img_ch2 = np.random.randn(nZ, nX, nRep) + 1j * np.random.randn(nZ, nX, nRep)
    
    # 计算
    dopu_result = calculate_spatial_dopu(img_ch1, img_ch2, params, use_gpu=True)
    
    print(f"DOPU结果形状: {dopu_result.shape}")
    print(f"DOPU范围: [{dopu_result.min():.4f}, {dopu_result.max():.4f}]")
