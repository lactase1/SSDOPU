"""
组合DOPU计算模块
计算分裂谱+空间组合DOPU
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
from calculate_split_spectrum_dopu import calculate_split_spectrum_dopu


def calculate_combined_dopu(
    Bd1: np.ndarray,
    Bd2: np.ndarray,
    params,
    SPL: int,
    nX: int,
    nr: int,
    nWin: int,
    windex: np.ndarray,
    winL: int,
    winG: np.ndarray,
    czrg: np.ndarray,
    use_gpu: bool = True
) -> np.ndarray:
    """
    计算分裂谱+空间组合DOPU
    
    算法原理:
        1. 先进行分裂谱DOPU计算，得到每个像素的分裂谱DOPU
        2. 然后对分裂谱DOPU结果进行空间滤波(3x3邻域平均)
        3. 得到最终的组合DOPU结果
    
    Args:
        Bd1: 通道1的OCT信号矩阵
        Bd2: 通道2的OCT信号矩阵
        params: 参数结构体
        SPL: 信号长度
        nX: X方向像素数
        nr: 重复次数
        nWin: 分裂谱窗口数
        windex: 窗口起始索引数组
        winL: 单个窗口长度
        winG: 高斯窗口函数
        czrg: Z方向裁剪范围
        use_gpu: 是否使用GPU加速
    
    Returns:
        dopu_combined: 组合DOPU结果 [Z像素×X像素]
    """
    # 检查是否需要进行组合DOPU计算
    if not params.dopu.do_combined:
        nZcrop = len(czrg)
        return np.ones((nZcrop, nX))
    
    # 步骤1: 计算分裂谱DOPU
    dopu_split, _ = calculate_split_spectrum_dopu(
        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg, use_gpu
    )
    
    # 步骤2: 对分裂谱DOPU进行空间滤波
    nZ, nX = dopu_split.shape
    
    # 选择计算后端
    if use_gpu and GPU_AVAILABLE:
        xp = cp
        dopu_split = to_gpu(dopu_split)
    else:
        xp = np
    
    dopu_combined = xp.zeros((nZ, nX))
    
    # 3点均值核
    kernel = xp.ones(3) / 3
    
    # 对每个深度层进行空间滤波
    for z in range(nZ):
        # 获取当前层的分裂谱DOPU
        dopu_layer = dopu_split[z, :]
        
        # 边界扩充处理
        dopu_padded = xp.pad(dopu_layer, 1, mode='edge')
        
        # 使用3点均值滤波器进行空间平均
        dopu_filtered = xp.convolve(dopu_padded, kernel, mode='valid')
        
        # 存储结果
        dopu_combined[z, :] = dopu_filtered
    
    # 转回CPU
    if use_gpu and GPU_AVAILABLE:
        dopu_combined = to_cpu(dopu_combined)
    
    return dopu_combined


if __name__ == "__main__":
    from config_params import config_params
    from gpu_utils import tukey_window
    
    print("测试组合DOPU计算...")
    
    # 创建测试参数
    params = config_params()
    params.dopu.do_combined = True
    params.dopu.do_ssdopu = True
    
    # 创建测试数据
    SPL = 4096
    nX = 500
    nr = 4
    nWin = 9
    Blength = 2048
    winL = int(2 * Blength / (nWin + 1))
    windex = np.arange(1, Blength + 1, winL // 2)[:nWin]
    winG = tukey_window(winL, 0.25)
    czrg = np.arange(1, 321)
    
    Bd1 = np.random.randn(Blength, nX, nr).astype(np.float64)
    Bd2 = np.random.randn(Blength, nX, nr).astype(np.float64)
    
    # 计算
    dopu_result = calculate_combined_dopu(
        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg, use_gpu=True
    )
    
    print(f"DOPU结果形状: {dopu_result.shape}")
    print(f"DOPU范围: [{dopu_result.min():.4f}, {dopu_result.max():.4f}]")
