"""
分裂谱DOPU计算模块
计算分裂谱偏振均匀度(Degree of Polarization Uniformity)
"""

import numpy as np
from typing import Tuple, Optional

try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None

from gpu_utils import get_array_module, to_gpu, to_cpu, fft, tukey_window
from cumulative_quv import cumulative_quv


def calculate_split_spectrum_dopu(
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
) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算分裂谱偏振均匀度(DOPU)
    
    算法原理:
        - 将OCT信号频谱分割成多个子频带
        - 对每个子频带分别计算偏振参数(Stokes参数)
        - 通过频谱平均提高信噪比和空间分辨率
        - 计算DOPU = sqrt(<S1/S0>² + <S2/S0>² + <S3/S0>²)
    
    Args:
        Bd1: 通道1的OCT信号矩阵 [信号长度×X像素×重复次数]
        Bd2: 通道2的OCT信号矩阵 [信号长度×X像素×重复次数]
        params: 参数结构体
        SPL: 信号长度 (FFT点数)
        nX: X方向像素数
        nr: 重复次数
        nWin: 分裂谱窗口数 (默认9)
        windex: 窗口起始索引数组
        winL: 单个窗口长度
        winG: 高斯窗口函数
        czrg: Z方向裁剪范围
        use_gpu: 是否使用GPU加速
    
    Returns:
        dopu_splitSpectrum: 分裂谱DOPU结果 [Z像素×X像素]
        dopu_ss: 各重复次数的分裂谱DOPU [Z像素×X像素×重复次数]
    """
    # 参数验证
    if Bd1.shape != Bd2.shape:
        raise ValueError("Bd1和Bd2的尺寸必须相同")
    
    signal_length, nX_check, nr_check = Bd1.shape
    if signal_length != SPL or nX_check != nX or nr_check != nr:
        raise ValueError("输入信号尺寸与参数不匹配")
    
    # 计算Z方向裁剪后的长度
    nZcrop = len(czrg)
    
    # 检查是否需要进行分裂谱计算
    if not params.dopu.do_ssdopu:
        dopu_splitSpectrum = np.ones((nZcrop, nX))
        dopu_ss = np.ones((nZcrop, nX, nr))
        return dopu_splitSpectrum, dopu_ss
    
    # 选择计算后端
    if use_gpu and GPU_AVAILABLE:
        xp = cp
        Bd1 = to_gpu(Bd1)
        Bd2 = to_gpu(Bd2)
        winG = to_gpu(winG)
    else:
        xp = np
    
    verbose = getattr(params.dopu, 'verbose', False)
    if verbose:
        print(f"  分裂谱DOPU: 窗口数={nWin}, 重复={nr}, X={nX}, 裁剪Z={nZcrop}")
    
    # 步骤1: 创建数组存储分裂谱复数FFT结果
    Bimg1 = xp.zeros((SPL, nX, nr, nWin), dtype=xp.complex128)
    Bimg2 = xp.zeros((SPL, nX, nr, nWin), dtype=xp.complex128)
    
    # 步骤2: 对每个重复和每个窗口进行FFT
    for iR in range(nr):
        for iL in range(nWin):
            # 提取当前窗口的数据片段并应用高斯窗口
            start_idx = int(windex[iL]) - 1  # 转为0索引
            end_idx = min(start_idx + winL, signal_length)
            
            # 确保窗口长度正确
            current_winL = end_idx - start_idx
            if current_winL != winL and len(winG) == winL:
                current_winG = winG[:current_winL]
            else:
                current_winG = winG
            
            # 应用窗口函数
            iBd1 = Bd1[start_idx:end_idx, :, iR] * current_winG[:, None]
            iBd2 = Bd2[start_idx:end_idx, :, iR] * current_winG[:, None]
            
            # 执行FFT并存储结果
            fft_result1 = xp.fft.fft(iBd1, n=SPL, axis=0)
            fft_result2 = xp.fft.fft(iBd2, n=SPL, axis=0)
            
            Bimg1[:, :, iR, iL] = fft_result1
            Bimg2[:, :, iR, iL] = fft_result2
    
    # 步骤3: 裁剪Z方向范围
    czrg_idx = czrg - 1  # 转为0索引
    IMGs_ch1 = Bimg1[czrg_idx, :, :, :]
    IMGs_ch2 = Bimg2[czrg_idx, :, :, :]
    
    # 步骤4: 计算每个窗口和重复的Stokes参数
    S0 = xp.zeros((nZcrop, nX, nr, nWin))
    S1 = xp.zeros((nZcrop, nX, nr, nWin))
    S2 = xp.zeros((nZcrop, nX, nr, nWin))
    S3 = xp.zeros((nZcrop, nX, nr, nWin))
    
    for iR in range(nr):
        for iL in range(nWin):
            s0, s1, s2, s3 = cumulative_quv(IMGs_ch1[:, :, iR, iL], IMGs_ch2[:, :, iR, iL])
            S0[:, :, iR, iL] = s0
            S1[:, :, iR, iL] = s1
            S2[:, :, iR, iL] = s2
            S3[:, :, iR, iL] = s3
    
    # 步骤5: 计算分裂谱DOPU
    # 避免除零
    S0_safe = xp.maximum(S0, xp.finfo(float).eps)
    
    # 对所有窗口进行平均，然后计算DOPU
    S1_avg = xp.mean(S1 / S0_safe, axis=3)  # [nZcrop×nX×nr]
    S2_avg = xp.mean(S2 / S0_safe, axis=3)
    S3_avg = xp.mean(S3 / S0_safe, axis=3)
    
    # 计算DOPU: sqrt(<S1/S0>² + <S2/S0>² + <S3/S0>²)
    dopu_ss = xp.sqrt(S1_avg**2 + S2_avg**2 + S3_avg**2)
    
    # 对所有重复进行平均，得到最终结果
    dopu_splitSpectrum = xp.mean(dopu_ss, axis=2)
    
    if verbose:
        print("  分裂谱DOPU计算完成")
    
    # 转回CPU
    if use_gpu and GPU_AVAILABLE:
        dopu_splitSpectrum = to_cpu(dopu_splitSpectrum)
        dopu_ss = to_cpu(dopu_ss)
    
    return dopu_splitSpectrum, dopu_ss


if __name__ == "__main__":
    from config_params import config_params
    
    print("测试分裂谱DOPU计算...")
    
    # 创建测试参数
    params = config_params()
    params.dopu.verbose = True
    
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
    dopu_result, dopu_all = calculate_split_spectrum_dopu(
        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg, use_gpu=True
    )
    
    print(f"DOPU结果形状: {dopu_result.shape}")
    print(f"DOPU范围: [{dopu_result.min():.4f}, {dopu_result.max():.4f}]")
