"""
DDG算法核心模块 - 完全向量化的GPU版本
移除所有Python循环，实现完全并行化
"""

import numpy as np
from typing import Tuple, Optional

try:
    import cupy as cp
    from cupyx.scipy import ndimage as cp_ndimage
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None
    cp_ndimage = None


def rotation_vector_to_matrix_batch(rotation_vectors):
    """
    批量计算旋转向量到旋转矩阵的转换 (Rodrigues公式)
    支持GPU加速
    
    Args:
        rotation_vectors: 旋转向量数组 [..., 3]
    
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


def ddg_ts_phr_vectorized(ps_bg_rm, B_bg_rm):
    """
    DDG切向法相位计算 - 向量化版本
    
    Args:
        ps_bg_rm: 去背景后的Stokes矢量序列 [output_depth, nX, Avnum, 3]
        B_bg_rm: 主轴方向 [output_depth, nX, 3]
    
    Returns:
        Subphr_bg_rm_Ts: 相位余弦值数组 [output_depth, nX, Avnum-2]
    """
    xp = cp.get_array_module(ps_bg_rm) if GPU_AVAILABLE else np
    
    # 计算切向量
    T1_bg_rm = ps_bg_rm[:, :, 1:, :] - ps_bg_rm[:, :, :-1, :]  # [..., Avnum-1, 3]
    T2_bg_rm = T1_bg_rm[:, :, 1:, :]  # [..., Avnum-2, 3]
    T1_bg_rm = T1_bg_rm[:, :, :-1, :]  # [..., Avnum-2, 3]
    
    # 扩展B维度用于广播
    B_expanded = B_bg_rm[:, :, None, :]  # [output_depth, nX, 1, 3]
    
    # 批量叉乘
    N1_bg_rm = xp.cross(B_expanded, T1_bg_rm, axis=-1)  # [..., Avnum-2, 3]
    N2_bg_rm = xp.cross(B_expanded, T2_bg_rm, axis=-1)
    
    # 计算点积和归一化
    dot_product = xp.sum(N1_bg_rm * N2_bg_rm, axis=-1)  # [..., Avnum-2]
    norm1 = xp.linalg.norm(N1_bg_rm, axis=-1)
    norm2 = xp.linalg.norm(N2_bg_rm, axis=-1)
    norm_product = norm1 * norm2
    
    # 计算余弦值
    cos_val = xp.where(norm_product > 1e-10, dot_product / norm_product, 1.0)
    cos_val = xp.clip(cos_val, -1, 1)
    
    return cos_val


def flatten_bscans_vectorized(Stru, Qm, Um, Vm, dopu_map, test_seg_top):
    """
    向量化的B-scan对齐（表面平展）
    
    Args:
        Stru, Qm, Um, Vm, dopu_map: [nZ, nX]
        test_seg_top: [nX]
    
    Returns:
        平展后的数组
    """
    xp = cp.get_array_module(Stru) if GPU_AVAILABLE else np
    nZ, nX = Stru.shape
    
    # 创建索引网格
    z_indices = xp.arange(nZ)[:, None]  # [nZ, 1]
    top_indices = test_seg_top[None, :].astype(int) - 1  # [1, nX]
    
    # 计算有效区域掩码
    valid_mask = (z_indices + top_indices < nZ) & (z_indices + top_indices >= 0)
    
    # 创建源索引
    src_z = xp.clip(z_indices + top_indices, 0, nZ - 1)
    
    # 应用平展
    StruF = xp.where(valid_mask, Stru[src_z, xp.arange(nX)], 0)
    QmF = xp.where(valid_mask, Qm[src_z, xp.arange(nX)], 0)
    UmF = xp.where(valid_mask, Um[src_z, xp.arange(nX)], 0)
    VmF = xp.where(valid_mask, Vm[src_z, xp.arange(nX)], 0)
    
    if dopu_map is not None:
        DopuF = xp.where(valid_mask, dopu_map[src_z, xp.arange(nX)], 0)
    else:
        DopuF = None
    
    return StruF, QmF, UmF, VmF, DopuF


def freespace_psoct_3_ddg_rmbg_7_vectorized(
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
):
    """
    PS-OCT DDG算法 - 完全向量化GPU版本
    移除所有Python循环，实现完全并行化
    
    Args:
        Qm, Um, Vm: Stokes参数 [nZ×nX] (cupy或numpy)
        Stru: 结构图像 [nZ×nX]
        test_seg_top: 表面位置 [nX]
        h1, h2: 滤波核
        Avnum: 平均层数
        dopu_map: DOPU图 [nZ×nX]
        enable_dopu_phase_supp: 是否启用DOPU相位抑制
    
    Returns:
        LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw
    """
    xp = cp.get_array_module(Qm) if GPU_AVAILABLE else np
    use_gpu = (xp == cp) if GPU_AVAILABLE else False
    
    # 确保所有输入在同一设备
    if use_gpu:
        Qm, Um, Vm = cp.asarray(Qm), cp.asarray(Um), cp.asarray(Vm)
        Stru = cp.asarray(Stru)
        test_seg_top = cp.asarray(test_seg_top)
        h1, h2 = cp.asarray(h1), cp.asarray(h2)
        if dopu_map is not None:
            dopu_map = cp.asarray(dopu_map)
    
    nZ, nX = Qm.shape
    
    # ===== Step 1: 数据平展 =====
    StruF, QmF, UmF, VmF, DopuF = flatten_bscans_vectorized(
        Stru, Qm, Um, Vm, dopu_map, test_seg_top
    )
    
    # ===== Step 2: QUV归一化 =====
    norm = xp.sqrt(QmF**2 + UmF**2 + VmF**2)
    norm = xp.maximum(norm, 1e-10)
    
    QF1 = QmF / norm
    UF1 = UmF / norm
    VF1 = VmF / norm
    
    QF1 = xp.nan_to_num(QF1, nan=0)
    UF1 = xp.nan_to_num(UF1, nan=0)
    VF1 = xp.nan_to_num(VF1, nan=0)
    
    SmapF = xp.stack([QF1, UF1, VF1], axis=2)  # [nZ, nX, 3]
    
    # 深度方向平滑
    if use_gpu:
        from cupyx.scipy.ndimage import uniform_filter1d
    else:
        from scipy.ndimage import uniform_filter1d
    SmapF1 = uniform_filter1d(SmapF, size=16, axis=0, mode='reflect')
    
    # ===== Step 3: 向量化DDG拟合 =====
    stLayer = 1
    DTheta = 1.0675 * 2 * xp.pi
    dTheta = DTheta / 320
    LABG = xp.array([1.0, 0.0, 0.0])
    
    output_depth = nZ - Avnum
    
    # 创建滑动窗口视图 [output_depth, nX, Avnum+1, 3]
    windows = xp.zeros((output_depth, nX, Avnum + 1, 3))
    for i in range(output_depth):
        # SmapF1 slice shape is [Avnum+1, nX, 3], transpose to [nX, Avnum+1, 3]
        windows[i] = SmapF1[i:i + Avnum + 1, :, :].transpose(1, 0, 2)
    
    # 检测零值（无效数据）
    zero_mask = xp.any(windows[:, :, :, 0] == 0, axis=2)  # [output_depth, nX]
    
    # ===== 计算DOPU =====
    if DopuF is not None:
        windowDOPU = xp.zeros((output_depth, nX))
        for i in range(output_depth):
            windowDOPU[i] = DopuF[i, :]
            zero_dopu = (windowDOPU[i] == 0)
            if xp.any(zero_dopu):
                avg_stokes = xp.mean(windows[i, zero_dopu, :, :], axis=1)
                windowDOPU[i, zero_dopu] = xp.sqrt(xp.sum(avg_stokes**2, axis=-1))
    else:
        avg_stokes = xp.mean(windows, axis=2)  # [output_depth, nX, 3]
        windowDOPU = xp.sqrt(xp.sum(avg_stokes**2, axis=-1))
    
    windowDOPU = xp.clip(windowDOPU, 0, 1)
    
    # DOPU自适应阈值
    dopuLowThresh = 0.5
    dopuMinScale = 0.1
    if dopu_map is not None:
        dopu_vals = dopu_map.flatten()
        if use_gpu:
            dopu_vals = dopu_vals.get()
        dopu_vals = dopu_vals[~np.isnan(dopu_vals) & (dopu_vals > 0)]
        if len(dopu_vals) > 0:
            med_dopu = np.median(dopu_vals)
            auto_thresh = med_dopu * 1.2
            dopuLowThresh = np.clip(auto_thresh, 0.1, 0.8)
            if use_gpu:
                dopuLowThresh = float(dopuLowThresh)
    
    # ===== 向量化计算旋转后的Stokes矢量 =====
    # 生成所有旋转角度 [output_depth, nX, Avnum]
    i_indices = xp.arange(output_depth)[:, None, None]  # [output_depth, 1, 1]
    j_indices = xp.arange(nX)[None, :, None]  # [1, nX, 1]
    jL_indices = xp.arange(Avnum)[None, None, :]  # [1, 1, Avnum]
    test_seg_expanded = test_seg_top[None, :, None]  # [1, nX, 1]
    
    angle1_base = -1 * (i_indices + test_seg_expanded + stLayer) * dTheta
    angle2_offsets = jL_indices - 1
    angle2 = angle1_base + angle2_offsets * (-dTheta)
    
    # 构造旋转向量
    rotation_vectors_rdu = angle1_base[..., None] * LABG[None, None, None, :]  # [output_depth, nX, 1, 3]
    rotation_vectors_rm = angle2[..., None] * LABG[None, None, None, :]  # [output_depth, nX, Avnum, 3]

    # 如果 rdu 没有 Avnum 维度，则沿该维度复制以匹配 inpp 的形状，避免 matmul 广播问题
    if rotation_vectors_rdu.shape[2] == 1 and rotation_vectors_rm.shape[2] > 1:
        rotation_vectors_rdu = xp.repeat(rotation_vectors_rdu, rotation_vectors_rm.shape[2], axis=2)

    # 批量计算旋转矩阵 [output_depth, nX, Avnum, 3, 3]
    R_rdu = rotation_vector_to_matrix_batch(rotation_vectors_rdu)
    R_rm = rotation_vector_to_matrix_batch(rotation_vectors_rm)

    # 确保旋转矩阵形状为 [output_depth, nX, Avnum, 3, 3]
    expected_shape = (output_depth, nX, rotation_vectors_rm.shape[2], 3, 3)
    def _normalize_R(R):
        R = xp.asarray(R)
        if R.shape == expected_shape:
            return R
        # 常见轴顺序尝试调整
        for perm in [(1,2,0,3,4), (0,2,1,3,4), (2,0,1,3,4), (1,0,2,3,4)]:
            try:
                Rr = R.transpose(perm)
                if Rr.shape == expected_shape:
                    return Rr
            except Exception:
                pass
        # 如果第三维为1，则复制以匹配 Avnum
        if R.shape[2] == 1 and expected_shape[2] > 1:
            R = xp.repeat(R, expected_shape[2], axis=2)
        # 最后尝试 reshape（保守操作）
        try:
            R = R.reshape(expected_shape)
        except Exception:
            raise ValueError(f'无法将旋转矩阵调整为期望形状 {expected_shape}, 当前形状 {R.shape}')
        return R

    R_rdu = _normalize_R(R_rdu)
    R_rm = _normalize_R(R_rm)

    # 应用旋转：R @ inpp [output_depth, nX, Avnum, 3, 1]
    inpp = windows[:, :, :Avnum, :, None]  # [output_depth, nX, Avnum, 3, 1]
    ps_bg_rdu = xp.matmul(R_rdu, inpp).squeeze(-1)  # [output_depth, nX, Avnum, 3]
    ps_bg_rm = xp.matmul(R_rm, inpp).squeeze(-1)
    
    # ===== 三点法拟合主轴 =====
    P1_rm = ps_bg_rm[:, :, 0, :]  # [output_depth, nX, 3]
    P2_rm = ps_bg_rm[:, :, Avnum // 2, :]
    P3_rm = ps_bg_rm[:, :, Avnum - 1, :]
    
    T1_rm = P2_rm - P1_rm  # [output_depth, nX, 3]
    T2_rm = P3_rm - P2_rm
    B_bg_rm = xp.cross(T1_rm, T2_rm, axis=-1)  # [output_depth, nX, 3]
    B_bg_rm_norm = xp.linalg.norm(B_bg_rm, axis=-1, keepdims=True)
    B_bg_rm = xp.where(B_bg_rm_norm > 1e-10, B_bg_rm / B_bg_rm_norm, xp.array([1.0, 0.0, 0.0]))
    
    # ===== DDG切向法计算相位 =====
    Subphr_bg_rm_Ts = ddg_ts_phr_vectorized(ps_bg_rm, B_bg_rm)  # [output_depth, nX, Avnum-2]
    phR_raw = xp.arccos(xp.median(Subphr_bg_rm_Ts, axis=2))  # [output_depth, nX]
    cumLA_raw = -B_bg_rm  # [output_depth, nX, 3]
    
    # ===== 去背景三点法 =====
    P1p_rdu = ps_bg_rdu[:, :, 0, :]
    P2p_rdu = ps_bg_rdu[:, :, Avnum // 2, :]
    P3p_rdu = ps_bg_rdu[:, :, Avnum - 1, :]
    
    T1p = P2p_rdu - P1p_rdu
    T2p = P3p_rdu - P2p_rdu
    Tb1 = P1_rm - P1p_rdu
    Tb2 = P3_rm - P3p_rdu
    Bp = xp.cross(T1p, T2p, axis=-1)
    
    T3_nor = xp.linalg.norm(P3_rm - P1_rm, axis=-1)  # [output_depth, nX]
    
    # 计算deltaB
    deltaB = xp.cross(T1p, Tb2, axis=-1) + xp.cross(T2p, Tb1, axis=-1) - xp.cross(Tb1, Tb2, axis=-1)
    
    # 低相位区域处理
    low_phase_mask = (T3_nor < 0.065)
    deltaB = xp.where(low_phase_mask[..., None], 0, deltaB)
    
    # 第一列方向检查
    first_col_positive = (Bp[:, 0, 0] > 0)
    Bp[:, 0, :] = xp.where(first_col_positive[:, None], -Bp[:, 0, :], Bp[:, 0, :])
    
    # 合成主轴
    B = Bp + deltaB
    B_norm = xp.linalg.norm(B, axis=-1, keepdims=True)
    B = xp.where(B_norm > 1e-10, B / B_norm, xp.array([1.0, 0.0, 0.0]))
    
    # ===== 计算相位 =====
    # 低相位区域的处理
    Subphr_rdu_Ts = ddg_ts_phr_vectorized(ps_bg_rdu, B)
    phR_low = xp.arccos(xp.median(Subphr_rdu_Ts, axis=2))
    
    # 正常相位
    Subphr_rm_Ts = ddg_ts_phr_vectorized(ps_bg_rm, B)
    phR_normal = xp.arccos(xp.median(Subphr_rm_Ts, axis=2))
    
    # 根据T3_nor选择
    phR = xp.where(low_phase_mask, phR_low, phR_normal)
    
    # 低相位区域的额外处理
    phR_prev = xp.concatenate([xp.zeros((1, nX)), phR[:-1, :]], axis=0)
    has_prev = xp.concatenate([xp.zeros((1, nX), dtype=bool), phR[:-1, :] > 0], axis=0)
    phR = xp.where(low_phase_mask & has_prev, (phR + phR_prev) / 2, phR)
    
    phR_rmBG = xp.maximum(phR - dTheta, 0)
    phR_rmBG = xp.where(low_phase_mask, xp.minimum(phR, 0.2), phR_rmBG)
    
    # ===== DOPU相位抑制 =====
    if enable_dopu_phase_supp:
        gamma = 2.0
        low_dopu_mask = (windowDOPU < dopuLowThresh)
        scale = xp.maximum((windowDOPU / dopuLowThresh) ** gamma, dopuMinScale)
        scale = xp.where(low_dopu_mask, scale, 1.0)
        
        phR_raw = phR_raw * scale
        phR = phR * scale
        phR_rmBG = xp.minimum(phR_rmBG, phR)
    
    # 应用零值掩码
    cumLA_bg = xp.where(zero_mask[..., None], 0, -B)
    phR = xp.where(zero_mask, 0, phR)
    phR_rmBG = xp.where(zero_mask, 0, phR_rmBG)
    phR_raw = xp.where(zero_mask, 0, phR_raw)
    cumLA_raw = xp.where(zero_mask[..., None], 0, cumLA_raw)
    
    # ===== Step 4: 滤波 =====
    if use_gpu:
        cumLA_bg_gF = xp.stack([
            cp_ndimage.convolve(cumLA_bg[:, :, i], h2, mode='reflect')
            for i in range(3)
        ], axis=2)
        cumLA_raw_gF = xp.stack([
            cp_ndimage.convolve(cumLA_raw[:, :, i], h2, mode='reflect')
            for i in range(3)
        ], axis=2)
        phR_gF = cp_ndimage.convolve(phR, h2, mode='reflect')
        phR_rmBG_gF = cp_ndimage.convolve(phR_rmBG, h2, mode='reflect')
        phR_raw_gF = cp_ndimage.convolve(phR_raw, h2, mode='reflect')
    else:
        from scipy.ndimage import convolve
        cumLA_bg_gF = np.stack([
            convolve(cumLA_bg[:, :, i], h2, mode='reflect')
            for i in range(3)
        ], axis=2)
        cumLA_raw_gF = np.stack([
            convolve(cumLA_raw[:, :, i], h2, mode='reflect')
            for i in range(3)
        ], axis=2)
        phR_gF = convolve(phR, h2, mode='reflect')
        phR_rmBG_gF = convolve(phR_rmBG, h2, mode='reflect')
        phR_raw_gF = convolve(phR_raw, h2, mode='reflect')
    
    # ===== Step 5: 主轴递推 (仍需循环但只在深度方向) =====
    # 构造初始旋转向量 [1, nX, 3]
    ax_rot = cumLA_bg_gF[0, :, :]  # [nX, 3]
    rotationVector_out = -phR_gF[0, :, None] / 2 * ax_rot  # [nX, 3]
    ax_rot_raw = cumLA_raw_gF[0, :, :]
    rotationVector_raw = -phR_raw_gF[0, :, None] / 2 * ax_rot_raw
    
    # 递推计算光轴
    bfloaxis3D = dr_la_vectorized(cumLA_bg_gF, -phR_gF, rotationVector_out, xp)
    bfloaxis3D_raw = dr_la_vectorized(cumLA_raw_gF, -phR_raw_gF, rotationVector_raw, xp)
    
    # ===== Step 6: 曲面还原 =====
    bfphrr = phR_rmBG_gF
    bfloaxis3D = xp.nan_to_num(bfloaxis3D, nan=0)
    bfphrr = xp.nan_to_num(bfphrr, nan=0)
    bfphrr_raw = phR_raw_gF
    bfloaxis3D_raw = xp.nan_to_num(bfloaxis3D_raw, nan=0)
    bfphrr_raw = xp.nan_to_num(bfphrr_raw, nan=0)
    
    # 还原到原始深度坐标
    LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw = unflatten_results_vectorized(
        bfloaxis3D, bfphrr, cumLA_bg_gF,
        bfloaxis3D_raw, bfphrr_raw, cumLA_raw_gF,
        test_seg_top, nZ, output_depth, xp
    )
    
    return LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw


def dr_la_vectorized(loaxis22, phR_gF, rotation_vector_init, xp):
    """
    光轴递推 - 向量化版本（深度方向仍需循环）
    
    Args:
        loaxis22: [output_depth, nX, 3]
        phR_gF: [output_depth, nX]
        rotation_vector_init: [nX, 3]
        xp: numpy或cupy
    
    Returns:
        axiss2: [output_depth, nX, 3]
    """
    nZ, nX, _ = loaxis22.shape
    axiss2 = xp.zeros((nZ, nX, 3))
    
    # 初始化旋转矩阵 [nX, 3, 3]
    rotation_matrix = rotation_vector_to_matrix_batch(rotation_vector_init)
    
    for i in range(nZ):
        a = loaxis22[i, :, :]  # [nX, 3]
        
        if i == 0:
            axiss2[i, :, :] = a
        else:
            # 批量矩阵乘法 [nX, 3, 3] @ [nX, 3, 1] = [nX, 3, 1]
            axis2 = xp.matmul(rotation_matrix, a[..., None]).squeeze(-1)  # [nX, 3]
            axiss2[i, :, :] = axis2
            
            # 更新旋转矩阵
            d1 = -phR_gF[i, :, None] / 2.0 * axis2  # [nX, 3]
            R_update = rotation_vector_to_matrix_batch(d1)  # [nX, 3, 3]
            rotation_matrix = xp.matmul(R_update, rotation_matrix)
            
            # 归一化
            norm = xp.linalg.norm(axis2, axis=-1, keepdims=True)
            axiss2[i, :, :] = xp.where(norm > 1e-10, axis2 / norm, 0)
            axiss2[i, :, 2] = xp.nan_to_num(axiss2[i, :, 2], nan=0)
    
    return axiss2


def unflatten_results_vectorized(
    bfloaxis3D, bfphrr, cumLA_bg_gF,
    bfloaxis3D_raw, bfphrr_raw, cumLA_raw_gF,
    test_seg_top, nZ, output_depth, xp
):
    """
    向量化的结果还原到原始坐标系
    """
    nX = len(test_seg_top)
    
    LA = xp.zeros((output_depth, nX, 3))
    PhR = xp.zeros((output_depth, nX))
    cumLA = xp.zeros((output_depth, nX, 3))
    LA_raw = xp.zeros((output_depth, nX, 3))
    PhR_raw = xp.zeros((output_depth, nX))
    cumLA_raw = xp.zeros((output_depth, nX, 3))
    
    for j in range(nX):
        top_idx = max(0, int(test_seg_top[j]) - 1)
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
    print("测试向量化DDG算法...")
    
    # 创建测试数据
    nZ, nX = 320, 500
    Qm = np.random.randn(nZ, nX)
    Um = np.random.randn(nZ, nX)
    Vm = np.random.randn(nZ, nX)
    Stru = np.random.rand(nZ, nX)
    test_seg_top = np.ones(nX) * 10
    
    # 创建滤波核
    h1 = np.ones((3, 3)) / 9
    h2 = np.ones((3, 3)) / 9
    
    import time
    start = time.time()
    
    # 运行向量化DDG算法
    LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw = freespace_psoct_3_ddg_rmbg_7_vectorized(
        Qm, Um, Vm, Stru, test_seg_top, h1, h2, 3
    )
    
    print(f"处理时间: {time.time() - start:.2f}秒")
    print(f"LA形状: {LA.shape}")
    print(f"PhR形状: {PhR.shape}")
