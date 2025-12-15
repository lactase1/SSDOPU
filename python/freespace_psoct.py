"""
DDG算法核心模块 - GPU向量化版本
PS-OCT后巩膜检测的核心算法实现
包括QUV归一化、DDG法光轴与相位拟合、光轴与相位滤波
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
    cp_ndimage = None

from gpu_utils import get_array_module, to_gpu, to_cpu, imfilter


def rotation_vector_to_matrix_batch(rotation_vectors):
    """
    批量计算旋转向量到旋转矩阵的转换 (Rodrigues公式)
    支持GPU加速
    
    Args:
        rotation_vectors: 旋转向量数组 [..., 3] (支持任意维度批量)
    
    Returns:
        rotation_matrices: 旋转矩阵数组 [..., 3, 3]
    """
    xp = cp.get_array_module(rotation_vectors) if GPU_AVAILABLE else np
    
    # 计算旋转角度
    theta = xp.linalg.norm(rotation_vectors, axis=-1, keepdims=True)  # [..., 1]
    
    # 处理小角度情况
    small_angle = (theta < 1e-10).squeeze(-1)
    
    # 归一化旋转轴
    k = xp.where(theta > 1e-10, rotation_vectors / theta, 0)  # [..., 3]
    
    # 构造反对称矩阵K [..., 3, 3]
    shape = rotation_vectors.shape[:-1]
    K = xp.zeros(shape + (3, 3), dtype=rotation_vectors.dtype)
    K[..., 0, 1] = -k[..., 2]
    K[..., 0, 2] = k[..., 1]
    K[..., 1, 0] = k[..., 2]
    K[..., 1, 2] = -k[..., 0]
    K[..., 2, 0] = -k[..., 1]
    K[..., 2, 1] = k[..., 0]
    
    # 计算 K @ K
    K2 = xp.matmul(K, K)
    
    # Rodrigues公式: R = I + sin(θ)K + (1-cos(θ))K²
    sin_theta = xp.sin(theta)[..., None]  # [..., 1, 1]
    cos_theta = xp.cos(theta)[..., None]  # [..., 1, 1]
    
    I = xp.eye(3, dtype=rotation_vectors.dtype)
    R = I + sin_theta * K + (1 - cos_theta) * K2
    
    # 小角度时返回单位矩阵
    R = xp.where(small_angle[..., None, None], I, R)
    
    return R


def rotation_vector_to_matrix(rotation_vector: np.ndarray) -> np.ndarray:
    """
    单个旋转向量转换为旋转矩阵 (保持向后兼容)
    
    Args:
        rotation_vector: 旋转向量 [3,]
    
    Returns:
        rotation_matrix: 3x3旋转矩阵
    """
    xp = cp.get_array_module(rotation_vector) if GPU_AVAILABLE else np
    result = rotation_vector_to_matrix_batch(rotation_vector[None, ...])
    return result[0]


def ddg_ts_phr(ps_bg_rm: np.ndarray, B_bg_rm: np.ndarray) -> np.ndarray:
    """
    DDG切向法相位计算
    
    Args:
        ps_bg_rm: 去背景后的Stokes矢量序列 [Avnum×3]
        B_bg_rm: 主轴方向 [3,]
    
    Returns:
        Subphr_bg_rm_Ts: 相位余弦值数组
    """
    Avnum = ps_bg_rm.shape[0]
    Subphr_bg_rm_Ts = np.zeros(Avnum - 2)
    
    for k in range(Avnum - 2):
        T1_bg_rm = ps_bg_rm[k + 1, :] - ps_bg_rm[k, :]
        T2_bg_rm = ps_bg_rm[k + 2, :] - ps_bg_rm[k + 1, :]
        N1_bg_rm = np.cross(B_bg_rm, T1_bg_rm)
        N2_bg_rm = np.cross(B_bg_rm, T2_bg_rm)
        
        norm_product = np.sqrt(np.sum(N1_bg_rm**2) * np.sum(N2_bg_rm**2))
        if norm_product > 1e-10:
            Subphr_bg_rm_Ts[k] = np.dot(N1_bg_rm, N2_bg_rm) / norm_product
        else:
            Subphr_bg_rm_Ts[k] = 1.0
    
    # 限制范围
    Subphr_bg_rm_Ts = np.clip(Subphr_bg_rm_Ts, -1, 1)
    return Subphr_bg_rm_Ts


def la_gf_filt(loaxis2: np.ndarray, h2: np.ndarray) -> np.ndarray:
    """
    光轴滤波（结构先验，空间一致性）
    
    Args:
        loaxis2: 光轴数据 [nZ×nX×3]
        h2: 滤波核
    
    Returns:
        loaxis3: 滤波后的光轴数据
    """
    Temp_ax1 = loaxis2[:, :, 0]
    Temp_ax2 = loaxis2[:, :, 1]
    Temp_ax3 = loaxis2[:, :, 2]
    
    Temp_ax1 = ndimage.convolve(Temp_ax1, h2, mode='reflect')
    Temp_ax2 = ndimage.convolve(Temp_ax2, h2, mode='reflect')
    Temp_ax3 = ndimage.convolve(Temp_ax3, h2, mode='reflect')
    
    # 归一化
    norm = np.sqrt(Temp_ax1**2 + Temp_ax2**2 + Temp_ax3**2)
    norm = np.maximum(norm, 1e-10)
    
    Temp_ax1 = Temp_ax1 / norm
    Temp_ax2 = Temp_ax2 / norm
    Temp_ax3 = Temp_ax3 / norm
    
    # 处理NaN
    Temp_ax1 = np.nan_to_num(Temp_ax1, nan=0)
    Temp_ax2 = np.nan_to_num(Temp_ax2, nan=0)
    Temp_ax3 = np.nan_to_num(Temp_ax3, nan=0)
    
    loaxis3 = np.stack([Temp_ax1, Temp_ax2, Temp_ax3], axis=2)
    return loaxis3


def dr_la(loaxis22: np.ndarray, phR_gF: np.ndarray, rotation_vector2: np.ndarray) -> np.ndarray:
    """
    光轴递推（主轴场恢复）
    
    Args:
        loaxis22: 光轴数据 [nZ×nX×3]
        phR_gF: 相位延迟 [nZ×nX]
        rotation_vector2: 初始旋转向量 [1×nX×3]
    
    Returns:
        axiss2: 递推后的光轴数据
    """
    nZ, nX, _ = loaxis22.shape
    axiss2 = np.zeros((nZ, nX, 3))
    axis2 = np.zeros((nZ, nX, 3))
    
    for j in range(nX):
        vector1 = rotation_vector2[0, j, :]
        rotation_matrix1 = rotation_vector_to_matrix(vector1)
        
        for i in range(nZ):
            a = loaxis22[i, j, :]
            
            if i == 0:
                axis2[i, j, :] = a
            else:
                axis2[i, j, :] = rotation_matrix1 @ a
                Qa, Ua, Va = axis2[i, j, :]
                d1 = -phR_gF[i, j] / 2.0 * np.array([Qa, Ua, Va])
                rotation_matrix1 = rotation_vector_to_matrix(d1) @ rotation_matrix1
                
                norm = np.sqrt(Qa**2 + Ua**2 + Va**2)
                if norm > 1e-10:
                    Q2 = Qa / norm
                    U2 = Ua / norm
                    V2 = Va / norm
                else:
                    Q2, U2, V2 = 0, 0, 0
                
                V2 = 0 if np.isnan(V2) else V2
                axiss2[i, j, :] = [Q2, U2, V2]
    
    return axiss2


def freespace_psoct_3_ddg_rmbg_7(
    Qm,
    Um,
    Vm,
    Stru,
    test_seg_top,
    h1,
    h2,
    Avnum: int,
    dopu_map=None,
    enable_dopu_phase_supp: bool = True
) -> Tuple:
    """
    PS-OCT DDG算法核心实现 - GPU向量化版本
    
    Args:
        Qm, Um, Vm: Stokes参数Q, U, V分量 [nZ×nX] (cupy或numpy数组)
        Stru: 结构图像 [nZ×nX]
        test_seg_top: 表面位置 [nX,]
        h1, h2: 滤波核
        Avnum: 平均层数
        dopu_map: DOPU图 [nZ×nX], 可选
        enable_dopu_phase_supp: 是否启用DOPU相位抑制
    
    Returns:
        LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw (全部为cupy或numpy数组)
    """
    # 判断使用GPU还是CPU
    xp = cp.get_array_module(Qm) if GPU_AVAILABLE else np
    use_gpu = (xp == cp) if GPU_AVAILABLE else False
    
    # 确保所有输入在同一设备
    if use_gpu:
        Qm = cp.asarray(Qm)
        Um = cp.asarray(Um)
        Vm = cp.asarray(Vm)
        Stru = cp.asarray(Stru)
        test_seg_top = cp.asarray(test_seg_top)
        if dopu_map is not None:
            dopu_map = cp.asarray(dopu_map)
    else:
        Qm = np.asarray(Qm)
        Um = np.asarray(Um)
        Vm = np.asarray(Vm)
        Stru = np.asarray(Stru)
        test_seg_top = np.asarray(test_seg_top)
        if dopu_map is not None:
            dopu_map = np.asarray(dopu_map)
    
    nZ, nX = Qm.shape
    
    # Step 1. 数据平展（A-scan表面对齐）
    test_seg_top = test_seg_top.astype(int)
    StruF = np.zeros_like(Stru)
    QmF = np.zeros_like(Qm)
    UmF = np.zeros_like(Um)
    VmF = np.zeros_like(Vm)
    DopuF = np.zeros_like(Stru) if dopu_map is not None else None
    
    for j in range(nX):
        top_idx = max(0, test_seg_top[j] - 1)  # 转为0索引
        valid_len = nZ - top_idx
        if valid_len > 0:
            StruF[:valid_len, j] = Stru[top_idx:, j]
            QmF[:valid_len, j] = Qm[top_idx:, j]
            UmF[:valid_len, j] = Um[top_idx:, j]
            VmF[:valid_len, j] = Vm[top_idx:, j]
            if dopu_map is not None:
                src_len = min(valid_len, dopu_map.shape[0] - top_idx)
                if src_len > 0:
                    DopuF[:src_len, j] = dopu_map[top_idx:top_idx + src_len, j]
    
    # Step 2. QUV归一化
    norm = np.sqrt(QmF**2 + UmF**2 + VmF**2)
    norm = np.maximum(norm, 1e-10)
    
    QF1 = QmF / norm
    UF1 = UmF / norm
    VF1 = VmF / norm
    
    QF1 = np.nan_to_num(QF1, nan=0)
    UF1 = np.nan_to_num(UF1, nan=0)
    VF1 = np.nan_to_num(VF1, nan=0)
    
    SmapF = np.stack([QF1, UF1, VF1], axis=2)
    
    # 深度方向平滑
    from scipy.ndimage import uniform_filter1d
    SmapF1 = uniform_filter1d(SmapF, size=16, axis=0, mode='reflect')
    
    # Step 3. DDG法光轴与相位拟合
    stLayer = 1
    DTheta = 1.0675 * 2 * np.pi
    dTheta = DTheta / 320
    LABG = np.array([1, 0, 0])
    
    output_depth = nZ - Avnum
    cumLA_bg = np.zeros((output_depth, nX, 3))
    phR = np.zeros((output_depth, nX))
    phR_rmBG = np.zeros((output_depth, nX))
    phR_raw_arr = np.zeros((output_depth, nX))
    cumLA_raw_arr = np.zeros((output_depth, nX, 3))
    phd = np.zeros((nZ, nX))
    rotationVector = np.zeros((output_depth, nX, 3))
    
    # DOPU阈值设置
    dopuLowThresh = 0.5
    dopuMinScale = 0.1
    
    if dopu_map is not None:
        dopu_vals = dopu_map.flatten()
        dopu_vals = dopu_vals[~np.isnan(dopu_vals) & (dopu_vals > 0)]
        if len(dopu_vals) > 0:
            med_dopu = np.median(dopu_vals)
            auto_ratio = 1.2
            auto_thresh = med_dopu * auto_ratio
            dopuLowThresh = np.clip(auto_thresh, 0.1, 0.8)
    
    # DDG拟合主循环
    for j in range(nX):
        dirCheck = True
        for i in range(output_depth):
            planeData = SmapF1[i:i + Avnum + 1, j, :]
            x1 = planeData[:, 0]
            
            if np.any(x1 == 0) or np.sum(x1 == 0) > 0:
                cumLA_bg[i, j, :] = 0
                phR[i, j] = 0
                phR_rmBG[i, j] = 0
                phR_raw_arr[i, j] = 0
                phd[i, j] = 0
            else:
                # 计算DOPU
                if DopuF is not None:
                    windowDOPU = DopuF[i, j]
                    if windowDOPU == 0:
                        avg_stokes = np.mean(planeData, axis=0)
                        windowDOPU = np.sqrt(np.sum(avg_stokes**2))
                else:
                    avg_stokes = np.mean(planeData, axis=0)
                    windowDOPU = np.sqrt(np.sum(avg_stokes**2))
                
                windowDOPU = np.clip(windowDOPU, 0, 1)
                phd[i, j] = windowDOPU
                
                # 计算旋转后的Stokes矢量
                ps_bg_rdu = []
                ps_bg_rm = []
                
                for jL in range(Avnum):
                    inpp = planeData[jL, :]
                    angle1 = -1 * (i + test_seg_top[j] + stLayer) * dTheta
                    angle2 = -1 * (i + test_seg_top[j] + stLayer + jL - 1) * dTheta
                    
                    outp = rotation_vector_to_matrix(angle1 * LABG) @ inpp
                    outp2 = rotation_vector_to_matrix(angle2 * LABG) @ inpp
                    
                    ps_bg_rdu.append(outp)
                    ps_bg_rm.append(outp2)
                
                ps_bg_rdu = np.array(ps_bg_rdu)
                ps_bg_rm = np.array(ps_bg_rm)
                
                # 三点法拟合主轴
                P1 = ps_bg_rm[0, :]
                P2 = ps_bg_rm[Avnum // 2, :]
                P3 = ps_bg_rm[Avnum - 1, :]
                
                T1 = P2 - P1
                T2 = P3 - P2
                B_bg_rm = np.cross(T1, T2)
                B_bg_rm_norm = np.linalg.norm(B_bg_rm)
                if B_bg_rm_norm > 1e-10:
                    B_bg_rm = B_bg_rm / B_bg_rm_norm
                else:
                    B_bg_rm = np.array([1, 0, 0])
                
                # DDG切向法计算相位
                Subphr_bg_rm_Ts = ddg_ts_phr(ps_bg_rm, B_bg_rm)
                phR_raw_arr[i, j] = np.arccos(np.median(Subphr_bg_rm_Ts))
                cumLA_raw_arr[i, j, :] = -B_bg_rm
                
                # 去背景三点法
                P1p = ps_bg_rdu[0, :]
                P2p = ps_bg_rdu[Avnum // 2, :]
                P3p = ps_bg_rdu[Avnum - 1, :]
                
                T1p = P2p - P1p
                T2p = P3p - P2p
                Tb1 = P1 - P1p
                Tb2 = P3 - P3p
                Bp = np.cross(T1p, T2p)
                T3_nor = np.linalg.norm(P3 - P1)
                
                if T3_nor < 0.065:
                    deltaB = np.array([0, 0, 0])
                    if dirCheck and Bp[0] > 0:
                        Bp = -Bp
                else:
                    dirCheck = False
                    deltaB = np.cross(T1p, Tb2) + np.cross(T2p, Tb1) - np.cross(Tb1, Tb2)
                
                B = Bp + deltaB
                B_norm = np.linalg.norm(B)
                if B_norm > 1e-10:
                    B = B / B_norm
                else:
                    B = np.array([1, 0, 0])
                
                if T3_nor < 0.065:
                    phR[i, j] = np.arccos(np.median(ddg_ts_phr(ps_bg_rdu, B)))
                    if i > 0 and phR[i - 1, j] > 0:
                        phR[i, j] = (phR[i, j] + phR[i - 1, j]) / 2
                    phR_rmBG[i, j] = max(phR[i, j] - dTheta, 0)
                    phR_rmBG[i, j] = min(phR[i, j], 0.2)
                else:
                    phR[i, j] = np.arccos(np.median(ddg_ts_phr(ps_bg_rm, B)))
                    phR_rmBG[i, j] = phR[i, j]
                
                # DOPU相位抑制
                if enable_dopu_phase_supp and windowDOPU < dopuLowThresh:
                    gamma = 2.0
                    scale = max((windowDOPU / dopuLowThresh) ** gamma, dopuMinScale)
                    phR_raw_arr[i, j] *= scale
                    phR[i, j] *= scale
                    phR_rmBG[i, j] = min(phR_rmBG[i, j], phR[i, j])
                
                cumLA_bg[i, j, :] = -B
    
    # Step 4. 光轴与相位滤波
    cumLA_bg_gF = la_gf_filt(cumLA_bg, h2)
    cumLA_raw_gF = la_gf_filt(cumLA_raw_arr, h2)
    phR_gF = ndimage.convolve(phR, h2, mode='reflect')
    phR_rmBG_gF = ndimage.convolve(phR_rmBG, h2, mode='reflect')
    phR_raw_gF = ndimage.convolve(phR_raw_arr, h2, mode='reflect')
    
    # Step 5. 旋转矩阵与主轴递推
    rotationVector_out = np.zeros((1, nX, 3))
    rotationVector_raw = np.zeros((1, nX, 3))
    
    for j in range(nX):
        ax_rot = cumLA_bg_gF[0, j, :]
        rotationVector_out[0, j, :] = -phR_gF[0, j] / 2 * ax_rot
        ax_rot_raw = cumLA_raw_gF[0, j, :]
        rotationVector_raw[0, j, :] = -phR_raw_gF[0, j] / 2 * ax_rot_raw
    
    bfloaxis3D = dr_la(cumLA_bg_gF, -phR_gF, rotationVector_out)
    bfloaxis3D_raw = dr_la(cumLA_raw_gF, -phR_raw_gF, rotationVector_raw)
    
    # Step 6. 曲面还原
    bfphrr = phR_rmBG_gF
    bfloaxis3D = np.nan_to_num(bfloaxis3D, nan=0)
    bfphrr = np.nan_to_num(bfphrr, nan=0)
    bfphrr_raw = phR_raw_gF
    bfloaxis3D_raw = np.nan_to_num(bfloaxis3D_raw, nan=0)
    bfphrr_raw = np.nan_to_num(bfphrr_raw, nan=0)
    
    LA = np.zeros((output_depth, nX, 3))
    PhR = np.zeros((output_depth, nX))
    cumLA = np.zeros((output_depth, nX, 3))
    LA_raw = np.zeros((output_depth, nX, 3))
    PhR_raw = np.zeros((output_depth, nX))
    cumLA_raw = np.zeros((output_depth, nX, 3))
    
    for j in range(nX):
        top_idx = max(0, test_seg_top[j] - 1)
        to_len = output_depth - top_idx
        from_len = min(to_len, output_depth)
        
        if from_len > 0:
            LA[top_idx:top_idx + from_len, j, :] = bfloaxis3D[:from_len, j, :]
            PhR[top_idx:top_idx + from_len, j] = bfphrr[:from_len, j]
            cumLA[top_idx:top_idx + from_len, j, :] = cumLA_bg_gF[:from_len, j, :]
            LA_raw[top_idx:top_idx + from_len, j, :] = bfloaxis3D_raw[:from_len, j, :]
            PhR_raw[top_idx:top_idx + from_len, j] = bfphrr_raw[:from_len, j]
            cumLA_raw[top_idx:top_idx + from_len, j, :] = cumLA_raw_gF[:from_len, j, :]
    
    return LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw


if __name__ == "__main__":
    print("测试DDG算法...")
    
    # 创建测试数据
    nZ, nX = 320, 500
    Qm = np.random.randn(nZ, nX)
    Um = np.random.randn(nZ, nX)
    Vm = np.random.randn(nZ, nX)
    Stru = np.random.rand(nZ, nX)
    test_seg_top = np.ones(nX) * 10
    
    # 创建滤波核
    from config_params import config_params
    params = config_params()
    h1 = params.filters.h1
    h2 = params.filters.h2
    
    # 运行DDG算法
    LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw = freespace_psoct_3_ddg_rmbg_7(
        Qm, Um, Vm, Stru, test_seg_top, h1, h2, params.polarization.avnum
    )
    
    print(f"LA形状: {LA.shape}")
    print(f"PhR形状: {PhR.shape}")
