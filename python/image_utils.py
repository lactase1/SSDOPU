"""
图像处理工具模块
包含表面分割、OAC计算、颜色映射等功能
"""

import numpy as np
from typing import Tuple, Optional
from scipy import ndimage
from scipy.signal import find_peaks, medfilt

try:
    import cupy as cp
    from cupyx.scipy import ndimage as cp_ndimage
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None


def cal_oac(lin_frame: np.ndarray) -> np.ndarray:
    """
    计算光学衰减系数(Optical Attenuation Coefficient)
    
    基于Beer-Lambert定律: I(z) = I₀ exp(-2μz)
    μ为衰减系数，反映组织的散射和吸收特性
    
    Args:
        lin_frame: 线性OCT信号强度矩阵 [Z×X]
    
    Returns:
        OAC: 光学衰减系数矩阵 [Z×X]
    """
    OAC = np.zeros_like(lin_frame)
    
    # 使用最后3层估算噪声底限
    last3 = lin_frame[-3:, :]
    Nois_M = np.mean(last3)
    
    if Nois_M < 1:
        return OAC
    
    # 噪声标准差
    Nois_D = np.std(last3)
    
    # 去除噪声底限并添加最小阈值
    lin_frame = lin_frame - Nois_M + Nois_D
    lin_frame = np.maximum(lin_frame, 1)  # 避免对数或分母为零
    
    # 估算尾部残余光强 (最后6层)
    tail_signal = np.mean(lin_frame[-6:, :], axis=0)
    tail_signal = medfilt(tail_signal, 15)
    
    # 平滑
    from scipy.ndimage import uniform_filter1d
    tail_signal = uniform_filter1d(tail_signal, 15, mode='reflect')
    
    # 经验衰减参数 alpha
    alpha = 0.0086
    
    # 计算每一深度层的OAC
    Z = lin_frame.shape[0]
    for z in range(Z):
        denom = 2 * alpha * np.sum(lin_frame[z+1:Z, :], axis=0) + tail_signal
        OAC[z, :] = lin_frame[z, :] / np.maximum(denom, 1e-10)
    
    return OAC


def surf_seg(OAC: np.ndarray, surf_threshold: float = 0.25) -> np.ndarray:
    """
    基于OAC数据进行组织表面检测和分割
    
    算法流程:
        1. 使用梯度模板检测表面边界
        2. 形态学处理去除噪声
        3. 峰值检测确定表面位置
        4. 中值滤波和平滑处理获得连续表面
    
    Args:
        OAC: 光学衰减系数矩阵 [Z×X]
        surf_threshold: 表面检测阈值
    
    Returns:
        test_seg_top: 每个A-line的组织表面位置向量 [X,]
    """
    nZ, nX = OAC.shape
    test_seg_top = np.ones(nX)
    
    if np.sum(OAC) < 1:
        return test_seg_top
    
    # 定义上表面检测的梯度模板
    h1 = np.array([-0.2, -0.2, -0.2, -0.2, -0.2, 1]).reshape(-1, 1)
    
    # 使用梯度模板进行卷积
    OAC_G = ndimage.convolve(OAC, h1, mode='reflect')
    OAC_G = (OAC_G > surf_threshold).astype(float)
    
    # 形态学处理
    from skimage.morphology import disk, binary_dilation, binary_erosion, remove_small_objects
    
    se = disk(2)
    OAC_G = binary_dilation(OAC_G.astype(bool), se)
    OAC_G = binary_erosion(OAC_G, se)
    
    # 去除小连通区域
    OAC_G = remove_small_objects(OAC_G.astype(bool), min_size=20)
    OAC_G = OAC_G.astype(float)
    
    # 逐列检测表面位置
    ILM = np.ones(nX)
    Temp_locs = 1
    
    for j in range(nX):
        # 检查前3行是否有信号
        if np.sum(OAC_G[:3, j]) > 1:
            ILM[j] = 1
            Temp_locs = 1
            continue
        
        # 使用峰值检测
        peaks, _ = find_peaks(OAC_G[:, j], height=0.2)
        
        if len(peaks) == 0:
            peaks = [Temp_locs]
        
        ILM[j] = peaks[0] + 1
        Temp_locs = peaks[0]
    
    # 中值滤波
    ILM = medfilt(ILM.astype(float), 11)
    
    # 平滑滤波
    from scipy.ndimage import uniform_filter1d
    test_seg_top = np.round(uniform_filter1d(ILM, 15, mode='reflect')) + 1
    
    return test_seg_top


def qu_coloring(fLAs: np.ndarray, cmapshift: int = 0) -> np.ndarray:
    """
    将双折射的Q、U分量转换为HSV彩色编码图像
    
    Args:
        fLAs: 双折射数据矩阵 [nZ×nX×3], 第3维为[?, Q, U]分量
        cmapshift: 颜色映射旋转角度
    
    Returns:
        hsvLA: HSV彩色编码图像 [nZ×nX×3]
    """
    import matplotlib.pyplot as plt
    
    # 创建HSV颜色映射表
    cmap = plt.cm.hsv(np.linspace(0, 1, 512))[:, :3]
    
    # 循环移位颜色表
    if cmapshift != 0:
        cmap = np.roll(cmap, cmapshift, axis=0)
    
    # 计算双折射轴角度
    thetas = np.arctan2(fLAs[:, :, 1], fLAs[:, :, 0])
    
    # 将角度映射到索引
    theta_range = np.linspace(-np.pi, np.pi, 512)
    theta_inds = np.interp(thetas.flatten(), theta_range, np.arange(512))
    theta_inds = np.round(theta_inds).astype(int)
    theta_inds = np.clip(theta_inds, 0, 511)
    
    # 获取颜色
    colors = cmap[theta_inds]
    hsvLA = colors.reshape((*fLAs.shape[:2], 3))
    
    return hsvLA


def mat2gray(arr: np.ndarray, limits: Optional[Tuple[float, float]] = None) -> np.ndarray:
    """
    将数组归一化到[0, 1]范围
    
    Args:
        arr: 输入数组
        limits: (min, max) 限制范围，默认使用数组的最小最大值
    
    Returns:
        归一化后的数组
    """
    if limits is None:
        limits = (arr.min(), arr.max())
    
    arr_out = (arr - limits[0]) / (limits[1] - limits[0] + 1e-10)
    return np.clip(arr_out, 0, 1)


def ind2rgb(ind: np.ndarray, cmap: np.ndarray) -> np.ndarray:
    """
    将索引图像转换为RGB图像
    
    Args:
        ind: 索引图像
        cmap: 颜色映射表 [N×3]
    
    Returns:
        RGB图像
    """
    ind = np.clip(ind.astype(int), 0, len(cmap) - 1)
    return cmap[ind]


def parula_colormap(n: int = 256) -> np.ndarray:
    """
    生成类似MATLAB parula的颜色映射
    
    Args:
        n: 颜色数量
    
    Returns:
        颜色映射表 [n×3]
    """
    # Parula颜色定义点
    parula_colors = np.array([
        [0.2422, 0.1504, 0.6603],
        [0.2810, 0.3228, 0.9579],
        [0.1786, 0.5289, 0.9682],
        [0.0689, 0.6948, 0.8394],
        [0.2161, 0.7843, 0.5923],
        [0.6720, 0.7793, 0.2227],
        [0.9970, 0.7659, 0.2199],
        [0.9769, 0.9839, 0.0805]
    ])
    
    # 插值生成完整的颜色映射
    from scipy.interpolate import interp1d
    
    x = np.linspace(0, 1, len(parula_colors))
    x_new = np.linspace(0, 1, n)
    
    cmap = np.zeros((n, 3))
    for i in range(3):
        f = interp1d(x, parula_colors[:, i], kind='cubic')
        cmap[:, i] = f(x_new)
    
    return np.clip(cmap, 0, 1)


if __name__ == "__main__":
    print("测试图像处理工具...")
    
    # 测试OAC计算
    nZ, nX = 320, 500
    lin_frame = np.random.rand(nZ, nX) * 1000
    oac = cal_oac(lin_frame)
    print(f"OAC形状: {oac.shape}")
    
    # 测试表面分割
    surf = surf_seg(oac)
    print(f"表面位置形状: {surf.shape}")
    
    # 测试颜色映射
    fLAs = np.random.randn(nZ, nX, 3)
    hsv_img = qu_coloring(fLAs)
    print(f"HSV图像形状: {hsv_img.shape}")
