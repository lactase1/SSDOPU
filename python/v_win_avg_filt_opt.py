"""
自适应DOPU滤波模块
根据DOPU值自适应调整高斯核大小进行滤波
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


def gaussian_kernel_2d(size: int, sigma: float) -> np.ndarray:
    """
    生成2D高斯核
    
    Args:
        size: 核大小
        sigma: 标准差
    
    Returns:
        高斯核
    """
    x = np.arange(size) - (size - 1) / 2
    y = np.arange(size) - (size - 1) / 2
    xx, yy = np.meshgrid(x, y)
    kernel = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
    return kernel / kernel.sum()


def v_win_avg_filt_opt(
    in_frame: np.ndarray,
    in_weight: np.ndarray,
    kRL: int,
    kRU: int,
    Nsec: Optional[int] = None,
    use_gpu: bool = True
) -> np.ndarray:
    """
    自适应DOPU滤波
    根据DOPU值(in_weight)自适应调整高斯核大小进行滤波
    
    Args:
        in_frame: 输入帧 [nZ×nX]
        in_weight: 置信度/DOPU值 [nZ×nX], 范围0~1
        kRL: 高斯核尺寸下限
        kRU: 高斯核尺寸上限
        Nsec: 核尺寸分段数
        use_gpu: 是否使用GPU加速
    
    Returns:
        out_frame: 滤波后的帧
    """
    if Nsec is None:
        Nsec = kRU - kRL + 1
    
    nZ, nX = in_frame.shape
    
    # 限制权重范围
    in_weight = np.clip(in_weight, 0, 1)
    in_weight = np.nan_to_num(in_weight, nan=0)
    
    if kRU <= 1:
        return in_frame * in_weight
    
    # 预生成核尺寸 (强度低 -> 大核，强度高 -> 小核)
    kRg = np.round(np.linspace(kRU, kRL, Nsec)).astype(int)
    kRg = kRg + (kRg + 1) % 2  # 调整为奇数便于对称窗口
    max_ker = kRg.max()
    max_rad = (max_ker - 1) // 2
    
    # 边界采用replicate padding
    in_frame_padded = np.pad(in_frame, max_rad, mode='edge')
    
    # 预生成所有需要的高斯核
    unique_sizes = np.unique(kRg)
    gauss_hs = {}
    for sz in unique_sizes:
        sigma = max(1, sz // 2)
        gauss_hs[sz] = gaussian_kernel_2d(sz, sigma)
    
    # 根据是否使用GPU选择计算方式
    if use_gpu and GPU_AVAILABLE:
        return _v_win_avg_filt_opt_gpu(
            in_frame, in_frame_padded, in_weight, kRg, Nsec, max_rad, gauss_hs
        )
    else:
        return _v_win_avg_filt_opt_cpu(
            in_frame, in_frame_padded, in_weight, kRg, Nsec, max_rad, gauss_hs
        )


def _v_win_avg_filt_opt_cpu(
    in_frame: np.ndarray,
    in_frame_padded: np.ndarray,
    in_weight: np.ndarray,
    kRg: np.ndarray,
    Nsec: int,
    max_rad: int,
    gauss_hs: dict
) -> np.ndarray:
    """CPU版本的自适应滤波"""
    nZ, nX = in_frame.shape
    out_frame = np.zeros_like(in_frame)
    
    for iz in range(nZ):
        for ix in range(nX):
            w_val = in_weight[iz, ix]
            # 计算使用哪个核尺寸
            f_ind = min(Nsec - 1, max(0, int(np.ceil(w_val * (Nsec - 1)))))
            ker_size = kRg[f_ind]
            rad = (ker_size - 1) // 2
            
            # 提取ROI
            rows = slice(iz + max_rad - rad, iz + max_rad + rad + 1)
            cols = slice(ix + max_rad - rad, ix + max_rad + rad + 1)
            roi = in_frame_padded[rows, cols]
            
            # 应用高斯核
            out_frame[iz, ix] = np.sum(roi * gauss_hs[ker_size])
    
    return out_frame


def _v_win_avg_filt_opt_gpu(
    in_frame: np.ndarray,
    in_frame_padded: np.ndarray,
    in_weight: np.ndarray,
    kRg: np.ndarray,
    Nsec: int,
    max_rad: int,
    gauss_hs: dict
) -> np.ndarray:
    """
    GPU加速版本的自适应滤波
    注意: 由于每个像素使用不同的核，无法完全并行化
    这里使用分组策略来加速
    """
    nZ, nX = in_frame.shape
    
    # 将数据转移到GPU
    in_frame_gpu = cp.asarray(in_frame)
    in_frame_padded_gpu = cp.asarray(in_frame_padded)
    in_weight_gpu = cp.asarray(in_weight)
    out_frame_gpu = cp.zeros_like(in_frame_gpu)
    
    # 将高斯核转移到GPU
    gauss_hs_gpu = {k: cp.asarray(v) for k, v in gauss_hs.items()}
    
    # 按核尺寸分组处理
    unique_sizes = np.unique(kRg)
    
    for ker_size in unique_sizes:
        rad = (ker_size - 1) // 2
        
        # 找到使用该核尺寸的像素索引
        f_ind_arr = cp.minimum(Nsec - 1, cp.maximum(0, cp.ceil(in_weight_gpu * (Nsec - 1)).astype(int)))
        mask = cp.zeros((nZ, nX), dtype=bool)
        for i, k in enumerate(kRg):
            if k == ker_size:
                mask |= (f_ind_arr == i)
        
        if not cp.any(mask):
            continue
        
        # 获取该核
        kernel = gauss_hs_gpu[ker_size]
        
        # 对这些像素进行滤波
        # 由于每个像素的核位置不同，需要逐个处理
        iz_arr, ix_arr = cp.where(mask)
        
        for idx in range(len(iz_arr)):
            iz = int(iz_arr[idx])
            ix = int(ix_arr[idx])
            
            rows = slice(iz + max_rad - rad, iz + max_rad + rad + 1)
            cols = slice(ix + max_rad - rad, ix + max_rad + rad + 1)
            roi = in_frame_padded_gpu[rows, cols]
            
            out_frame_gpu[iz, ix] = cp.sum(roi * kernel)
    
    return cp.asnumpy(out_frame_gpu)


def v_win_avg_filt_opt_vectorized(
    in_frame,
    in_weight,
    kRL: int,
    kRU: int,
    use_gpu: bool = True
):
    """
    向量化版本的自适应DOPU滤波 - GPU优化
    使用多个固定核的加权组合来近似自适应滤波
    
    Args:
        in_frame: 输入帧 [nZ×nX] (cupy或numpy数组)
        in_weight: 置信度/DOPU值 [nZ×nX], 范围0~1
        kRL: 高斯核尺寸下限
        kRU: 高斯核尺寸上限
        use_gpu: 是否使用GPU加速
    
    Returns:
        out_frame: 滤波后的帧 (与输入同类型)
    """
    # 判断使用GPU还是CPU
    use_gpu = use_gpu and GPU_AVAILABLE
    xp = cp.get_array_module(in_frame) if GPU_AVAILABLE else np
    
    # 确保数据在正确的设备上
    if use_gpu:
        in_frame = cp.asarray(in_frame)
        in_weight = cp.asarray(in_weight)
    else:
        in_frame = np.asarray(in_frame)
        in_weight = np.asarray(in_weight)
    
    # 限制权重范围
    in_weight = xp.clip(in_weight, 0, 1)
    in_weight = xp.nan_to_num(in_weight, nan=0)
    
    # 生成几个不同大小的滤波结果
    sigmas = xp.linspace(kRL / 2, kRU / 2, 5)
    filtered_results = []
    
    for sigma in sigmas:
        if use_gpu:
            from cupyx.scipy import ndimage as cp_ndimage
            filtered = cp_ndimage.gaussian_filter(in_frame, float(sigma))
        else:
            from scipy import ndimage
            filtered = ndimage.gaussian_filter(in_frame, float(sigma))
        filtered_results.append(filtered)
    
    # 根据权重插值
    # 权重高(接近1) -> 使用小核(index 0)
    # 权重低(接近0) -> 使用大核(index -1)
    out_frame = xp.zeros_like(in_frame)
    n_levels = len(filtered_results)
    
    for i in range(n_levels):
        # 计算每个level的权重
        level_weight = xp.clip(1 - xp.abs(in_weight - i / (n_levels - 1)) * (n_levels - 1), 0, 1)
        out_frame += level_weight * filtered_results[n_levels - 1 - i]
    
    # 归一化
    weight_sum = xp.zeros_like(in_frame)
    for i in range(n_levels):
        level_weight = xp.clip(1 - xp.abs(in_weight - i / (n_levels - 1)) * (n_levels - 1), 0, 1)
        weight_sum += level_weight
    
    out_frame = out_frame / xp.maximum(weight_sum, 1e-10)
    
    return out_frame


if __name__ == "__main__":
    print("测试自适应DOPU滤波...")
    
    # 创建测试数据
    nZ, nX = 320, 500
    in_frame = np.random.randn(nZ, nX)
    in_weight = np.random.rand(nZ, nX)
    
    # 测试CPU版本
    import time
    
    start = time.time()
    out_cpu = v_win_avg_filt_opt(in_frame, in_weight, 3, 21, use_gpu=False)
    print(f"CPU版本耗时: {time.time() - start:.2f}秒")
    
    # 测试向量化版本
    start = time.time()
    out_vec = v_win_avg_filt_opt_vectorized(in_frame, in_weight, 3, 21, use_gpu=False)
    print(f"向量化版本耗时: {time.time() - start:.2f}秒")
    
    # 测试GPU版本
    if GPU_AVAILABLE:
        start = time.time()
        out_gpu = v_win_avg_filt_opt(in_frame, in_weight, 3, 21, use_gpu=True)
        print(f"GPU版本耗时: {time.time() - start:.2f}秒")
