"""
PS-OCT主处理模块
整合所有处理步骤的主函数
"""

import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict, Any
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import time
import warnings

try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None

from config_params import config_params, ConfigParams
from oct_file_reader import OCTFileReader
from gpu_utils import to_gpu, to_cpu, tukey_window, hilbert
from cumulative_quv import cumulative_quv, cumulative_quv_gpu
from calculate_split_spectrum_dopu import calculate_split_spectrum_dopu
from calculate_spatial_dopu import calculate_spatial_dopu
from calculate_combined_dopu import calculate_combined_dopu
from v_win_avg_filt_opt import v_win_avg_filt_opt_vectorized
from freespace_psoct_vectorized import freespace_psoct_3_ddg_rmbg_7_vectorized
from image_utils import cal_oac, surf_seg, qu_coloring, mat2gray, parula_colormap
from dicom_utils import dicom_write, convert_dcm_to_tiff


def cal_la_phr_all(
    img_ch1,
    img_ch2,
    test_seg_top,
    dopu_split_spec_m,
    kRL: int,
    kRU: int,
    h1,
    h2,
    avnum: int,
    wov_win_f: int = 0,
    enable_dopu_phase_supp: bool = True
) -> Tuple:
    """
    计算完整的偏振参数 - GPU优化版本
    
    Args:
        img_ch1, img_ch2: 双通道复数OCT信号 (cupy或numpy数组)
        test_seg_top: 组织表面位置
        dopu_split_spec_m: 分光谱DOPU矩阵
        kRL, kRU: 滤波核范围参数
        h1, h2: 高斯滤波核
        avnum: 平均层数
        wov_win_f: 滤波模式 (1:固定高斯滤波, 0:自适应DOPU滤波)
        enable_dopu_phase_supp: 是否启用DOPU相位抑制
    
    Returns:
        LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw
    """
    # 判断是否在GPU上
    xp = cp.get_array_module(img_ch1) if GPU_AVAILABLE else np
    use_gpu = (xp == cp) if GPU_AVAILABLE else False
    
    # 确保数据在GPU上（如果可用）
    if use_gpu:
        img_ch1 = cp.asarray(img_ch1)
        img_ch2 = cp.asarray(img_ch2)
        test_seg_top = cp.asarray(test_seg_top)
        dopu_split_spec_m = cp.asarray(dopu_split_spec_m)
        h1, h2 = cp.asarray(h1), cp.asarray(h2)
    
    nZ, nX = img_ch1.shape
    
    # 计算Stokes参数（保持在GPU）
    if use_gpu:
        ES0, ES1, ES2, ES3 = cumulative_quv_gpu(img_ch1, img_ch2)
    else:
        ES0, ES1, ES2, ES3 = cumulative_quv(img_ch1, img_ch2)
    
    if xp.sum(ES0) < 5:
        # 返回零数组
        zero_shape = (nZ - avnum, nX, 3)
        return (xp.zeros(zero_shape), xp.zeros((nZ - avnum, nX)), xp.zeros(zero_shape),
                xp.zeros(zero_shape), xp.zeros((nZ - avnum, nX)), xp.zeros(zero_shape))
    
    # 计算归一化Stokes分量
    eps = float(xp.finfo(xp.float64).eps)
    ES0_safe = xp.maximum(ES0, eps)
    EQm = ES1 / ES0_safe
    EUm = ES2 / ES0_safe
    EVm = ES3 / ES0_safe
    
    # 初始化结构矩阵
    Stru_E = xp.zeros((nZ, nX))
    
    # 根据滤波模式进行预处理
    if wov_win_f == 1:
        # 模式1: 固定高斯滤波（在GPU上）
        if use_gpu:
            from cupyx.scipy import ndimage as cp_ndimage
            EQmm = cp_ndimage.convolve(EQm, h1, mode='reflect')
            EUmm = cp_ndimage.convolve(EUm, h1, mode='reflect')
            EVmm = cp_ndimage.convolve(EVm, h1, mode='reflect')
        else:
            from scipy import ndimage
            EQmm = ndimage.convolve(EQm, h1, mode='reflect')
            EUmm = ndimage.convolve(EUm, h1, mode='reflect')
            EVmm = ndimage.convolve(EVm, h1, mode='reflect')
    else:
        # 模式0: 自适应DOPU滤波（保持在GPU）
        EQmm = v_win_avg_filt_opt_vectorized(EQm, dopu_split_spec_m, kRL, kRU, use_gpu=use_gpu)
        EUmm = v_win_avg_filt_opt_vectorized(EUm, dopu_split_spec_m, kRL, kRU, use_gpu=use_gpu)
        EVmm = v_win_avg_filt_opt_vectorized(EVm, dopu_split_spec_m, kRL, kRU, use_gpu=use_gpu)
    
    # 调用向量化核心算法（数据保持在GPU）
    LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw = freespace_psoct_3_ddg_rmbg_7_vectorized(
        EQmm, EUmm, EVmm, Stru_E, test_seg_top, h1, h2, avnum, 
        dopu_split_spec_m, enable_dopu_phase_supp
    )
    
    return LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw


def process_bscan(
    iY: int,
    reader: OCTFileReader,
    params: ConfigParams,
    Ref_ch1: np.ndarray,
    Ref_ch2: np.ndarray,
    phV: np.ndarray,
    winG: np.ndarray,
    winG_whole: np.ndarray,
    nWin: int,
    windex: np.ndarray,
    winL: int,
    czrg: np.ndarray,
    topLines: np.ndarray,
    SPL: int,
    nX: int,
    nr: int,
    Blength: int
) -> Dict[str, Any]:
    """
    处理单个B-scan - GPU优化版本
    
    Args:
        iY: B-scan索引
        reader: OCT文件读取器
        params: 配置参数
        其他参数: 处理所需的各种参数
    
    Returns:
        包含处理结果的字典
    """
    # 读取B-scan数据
    Bs1, Bs2 = reader.read_bscan(iY, nr)
    
    # ===== 立即上传到GPU（如果可用）=====
    if GPU_AVAILABLE:
        Bs1_gpu = cp.asarray(Bs1)
        Bs2_gpu = cp.asarray(Bs2)
        Ref_ch1_gpu = cp.asarray(Ref_ch1[:Blength, :, np.newaxis])
        Ref_ch2_gpu = cp.asarray(Ref_ch2[:Blength, :, np.newaxis])
        phV_gpu = cp.asarray(phV[:, np.newaxis, np.newaxis])
        # winG和winG_whole保持一维，不添加额外维度
        winG_gpu = cp.asarray(winG)
        winG_whole_gpu = cp.asarray(winG_whole)
        czrg_gpu = cp.asarray(czrg)
        
        # 减去参考（在GPU上）
        Bs1_gpu = Bs1_gpu - Ref_ch1_gpu
        Bs2_gpu = Bs2_gpu - Ref_ch2_gpu
        
        # 色散补偿（在GPU上）
        if params.processing.do_ph_comp:
            from gpu_utils import hilbert as gpu_hilbert
            Bd1_gpu = cp.real(gpu_hilbert(Bs1_gpu, axis=0) * phV_gpu)
            Bd2_gpu = cp.real(gpu_hilbert(Bs2_gpu, axis=0) * phV_gpu)
        else:
            Bd1_gpu = Bs1_gpu
            Bd2_gpu = Bs2_gpu
    else:
        # CPU版本
        Bs1 = Bs1 - Ref_ch1[:Blength, :, np.newaxis]
        Bs2 = Bs2 - Ref_ch2[:Blength, :, np.newaxis]
        
        if params.processing.do_ph_comp:
            Bd1 = np.real(hilbert(Bs1, axis=0) * phV[:, np.newaxis, np.newaxis])
            Bd2 = np.real(hilbert(Bs2, axis=0) * phV[:, np.newaxis, np.newaxis])
        else:
            Bd1 = Bs1
            Bd2 = Bs2
        Bd1_gpu, Bd2_gpu = Bd1, Bd2
        winG_gpu, winG_whole_gpu = winG, winG_whole
        czrg_gpu = czrg
    
    nZcrop = len(czrg)
    
    # 计算DOPU（保持在GPU）
    if params.dopu.do_combined:
        dopu_result = calculate_combined_dopu(
            Bd1_gpu, Bd2_gpu, params, SPL, nX, nr, nWin, windex, winL, winG_gpu, czrg_gpu, 
            use_gpu=GPU_AVAILABLE
        )
        if GPU_AVAILABLE:
            dopu_ss = dopu_result[:, :, cp.newaxis]
        else:
            dopu_ss = dopu_result[:, :, np.newaxis]
    elif params.dopu.do_spatial:
        # 先计算FFT（在GPU上）
        if GPU_AVAILABLE:
            Bimg1_whole = cp.fft.fft(Bd1_gpu * winG_whole_gpu[:, cp.newaxis, cp.newaxis], n=SPL, axis=0)
            Bimg2_whole = cp.fft.fft(Bd2_gpu * winG_whole_gpu[:, cp.newaxis, cp.newaxis], n=SPL, axis=0)
            IMG1_whole = Bimg1_whole[czrg_gpu - 1, :, :]
            IMG2_whole = Bimg2_whole[czrg_gpu - 1, :, :]
        else:
            Bimg1_whole = np.fft.fft(Bd1_gpu * winG_whole_gpu[:, np.newaxis, np.newaxis], n=SPL, axis=0)
            Bimg2_whole = np.fft.fft(Bd2_gpu * winG_whole_gpu[:, np.newaxis, np.newaxis], n=SPL, axis=0)
            IMG1_whole = Bimg1_whole[czrg_gpu - 1, :, :]
            IMG2_whole = Bimg2_whole[czrg_gpu - 1, :, :]
        dopu_result = calculate_spatial_dopu(IMG1_whole, IMG2_whole, params, use_gpu=GPU_AVAILABLE)
        if GPU_AVAILABLE:
            dopu_ss = dopu_result[:, :, cp.newaxis]
        else:
            dopu_ss = dopu_result[:, :, np.newaxis]
    else:
        dopu_result, dopu_ss = calculate_split_spectrum_dopu(
            Bd1_gpu, Bd2_gpu, params, SPL, nX, nr, nWin, windex, winL, winG_gpu, czrg_gpu,
            use_gpu=GPU_AVAILABLE
        )
    
    # 全频谱FFT（在GPU上）
    if GPU_AVAILABLE:
        Bimg1_whole = cp.fft.fft(Bd1_gpu * winG_whole_gpu[:, cp.newaxis, cp.newaxis], n=SPL, axis=0)
        Bimg2_whole = cp.fft.fft(Bd2_gpu * winG_whole_gpu[:, cp.newaxis, cp.newaxis], n=SPL, axis=0)
        IMG1_whole = Bimg1_whole[czrg_gpu - 1, :, :]
        IMG2_whole = Bimg2_whole[czrg_gpu - 1, :, :]
        
        # 计算Stokes参数（在GPU上）
        wS0, wS1, wS2, wS3 = cumulative_quv_gpu(IMG1_whole, IMG2_whole)
    else:
        Bimg1_whole = np.fft.fft(Bd1_gpu * winG_whole_gpu[:, np.newaxis, np.newaxis], n=SPL, axis=0)
        Bimg2_whole = np.fft.fft(Bd2_gpu * winG_whole_gpu[:, np.newaxis, np.newaxis], n=SPL, axis=0)
        IMG1_whole = Bimg1_whole[czrg_gpu - 1, :, :]
        IMG2_whole = Bimg2_whole[czrg_gpu - 1, :, :]
        
        # 计算Stokes参数（在CPU上）
        wS0, wS1, wS2, wS3 = cumulative_quv(IMG1_whole, IMG2_whole)
    
    xp = cp if GPU_AVAILABLE else np
    eps = float(xp.finfo(xp.float64).eps)
    wS0_safe = xp.maximum(wS0, eps)
    wQ = wS1 / wS0_safe
    wU = wS2 / wS0_safe
    wV = wS3 / wS0_safe
    
    strLin = xp.mean(wS0, axis=2)
    Strus = 20 * xp.log10(xp.maximum(strLin, 1e-10))
    Smap_avg = xp.stack([xp.mean(wQ, axis=2), xp.mean(wU, axis=2), xp.mean(wV, axis=2)], axis=2)
    Smap_rep1 = xp.stack([wQ[:, :, 0], wU[:, :, 0], wV[:, :, 0]], axis=2)
    
    # 表面分割（需要在CPU上）
    if not params.processing.has_seg:
        if GPU_AVAILABLE:
            strLin_cpu = cp.asnumpy(strLin)
        else:
            strLin_cpu = strLin
        strOAC = cal_oac(strLin_cpu)
        topLines[:, iY] = surf_seg(strOAC, 0.25) + 2
    
    # 计算局部双折射和相位延迟（保持在GPU）
    dopu_split_spec_m = xp.mean(dopu_ss, axis=2)
    IMG_ch1 = xp.mean(IMG1_whole, axis=2)
    IMG_ch2 = xp.mean(IMG2_whole, axis=2)
    
    # 调用DDG算法（数据保持在GPU）
    if GPU_AVAILABLE:
        topLines_gpu = cp.asarray(topLines[:, iY])
    else:
        topLines_gpu = topLines[:, iY]
    
    LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw = cal_la_phr_all(
        IMG_ch1, IMG_ch2, topLines_gpu, dopu_split_spec_m,
        params.polarization.krl_cfg1, params.polarization.kru_cfg1,
        params.filters.h1, params.filters.h2, params.polarization.avnum,
        params.mode.wov_win_f, params.polarization.enable_dopu_phase_supp
    )
    
    # ===== 仅在最后转回CPU =====
    if GPU_AVAILABLE:
        dopu_result = cp.asnumpy(dopu_result)
        Strus = cp.asnumpy(Strus)
        Smap_avg = cp.asnumpy(Smap_avg)
        Smap_rep1 = cp.asnumpy(Smap_rep1)
        LA = cp.asnumpy(LA)
        PhR = cp.asnumpy(PhR)
        cumLA = cp.asnumpy(cumLA)
        LA_raw = cp.asnumpy(LA_raw)
        PhR_raw = cp.asnumpy(PhR_raw)
        cumLA_raw = cp.asnumpy(cumLA_raw)
    
    return {
        'iY': iY,
        'dopu': dopu_result,
        'Strus': Strus,
        'Smap_avg': Smap_avg,
        'Smap_rep1': Smap_rep1,
        'topLines': topLines[:, iY],
        'LA': LA,
        'PhR': PhR,
        'cumLA': cumLA,
        'LA_raw': LA_raw,
        'PhR_raw': PhR_raw,
        'cumLA_raw': cumLA_raw
    }


def process_single_file(input_file_path: str, output_base: str, params: Optional[ConfigParams] = None):
    """
    处理单个OCT文件
    
    Args:
        input_file_path: 输入文件路径
        output_base: 输出目录
        params: 配置参数，如果为None则使用默认参数
    """
    file_start_time = time.time()
    
    if params is None:
        params = config_params()
    
    input_path = Path(input_file_path)
    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file_path}")
    
    # 创建输出目录
    output_path = Path(output_base)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*50}")
    print(f"正在处理文件: {input_path.name}")
    print(f"文件路径: {input_file_path}")
    print(f"{'='*50}")
    
    # 读取文件头
    with OCTFileReader(input_file_path) as reader:
        header = reader.read_header()
        
        SPL = int(header['SPL'])
        nX = header['nX']
        nY = header['nY_actual']
        Blength = header['Blength']
        nR = header['nR']
        nr = nR
        disp_coef = header['disp_coef']
        
        print(f"SPL: {SPL}, nX: {nX}, nY: {nY}, Blength: {Blength}")
        
        # 帧数限制
        if params.processing.max_frames > 0:
            original_nY = nY
            nY = min(nY, params.processing.max_frames)
            print(f"帧数限制: {original_nY} -> {nY}")
        
        # 计算色散补偿相位
        kmat = np.linspace(-0.5, 0.5, Blength) ** 2
        phV = np.exp(1j * kmat * disp_coef)
        
        # 创建窗口函数
        nWin = 9
        winL = int(2 * Blength / (nWin + 1))
        winG = tukey_window(winL, 0.25)
        winG_whole = tukey_window(Blength, 0.25)
        windex = np.arange(1, Blength + 1, winL // 2)[:nWin]
        
        # 读取参考
        if params.processing.useref == 1:
            Ref_ch1, Ref_ch2 = reader.read_reference()
        elif params.processing.useref == -1:
            Ref_ch1 = np.zeros((SPL, nX))
            Ref_ch2 = np.zeros((SPL, nX))
        else:
            Ref_ch1 = np.tile(header['Bck1'][:, np.newaxis], (1, nX))
            Ref_ch2 = np.tile(header['Bck2'][:, np.newaxis], (1, nX))
        
        # 设置Z范围
        czrg = np.arange(1, 321)
        if params.range.set_z_rg > 0:
            czrg = czrg[:params.range.set_z_rg]
        nZcrop = len(czrg)
        
        # 初始化表面位置
        topLines = np.ones((nX, nY))
        
        # 加载预分割结果
        if params.processing.has_seg:
            mat_file = input_path.with_suffix('.mat')
            if mat_file.exists():
                try:
                    from scipy.io import loadmat
                    mat_data = loadmat(str(mat_file))
                    if 'topLines' in mat_data:
                        topLines = mat_data['topLines'].astype(float)
                        topLines[topLines <= 1] = 1
                        print(f"成功加载分割结果: {mat_file}")
                except Exception as e:
                    warnings.warn(f"加载分割结果失败: {e}")
                    params.processing.has_seg = False
            else:
                warnings.warn(f"分割结果文件不存在: {mat_file}")
                params.processing.has_seg = False
        
        # 创建输出目录
        name = input_path.stem
        if params.range.set_z_rg:
            foutputdir = output_path / f"{name}.rep{nr}fastZ_testDDG{params.range.set_z_rg}"
        else:
            foutputdir = output_path / f"{name}.rep{nr}DDG_corr3_wobg_test"
        foutputdir.mkdir(parents=True, exist_ok=True)
        
        # 初始化结果数组
        avnum = params.polarization.avnum
        LA_c_cfg1_avg = np.zeros((nZcrop - avnum, nX, 3, nY))
        PhR_c_cfg1_avg = np.zeros((nZcrop - avnum, nX, nY))
        cumLA_cfg1_avg = np.zeros((nZcrop - avnum, nX, 3, nY))
        LA_Ms_cfg1_rmBG = np.zeros_like(LA_c_cfg1_avg)
        PhR_Ms_cfg1_rmBG = np.zeros_like(PhR_c_cfg1_avg)
        cumLA_Ms_cfg1_rmBG = np.zeros_like(cumLA_cfg1_avg)
        Smap_avg = np.zeros((nZcrop, nX, 3, nY))
        Smap_rep1 = np.zeros((nZcrop, nX, 3, nY))
        Strus = np.zeros((nZcrop, nX, nY))
        dopu_splitSpectrum = np.zeros((nZcrop, nX, nY))
        
        print(f"开始处理 {nY} 个B-Scan...")
        
        # 处理每个B-scan
        for iY in range(nY):
            if (iY + 1) % 10 == 0 or iY == 0:
                print(f"  处理 B-scan {iY + 1}/{nY}...")
            
            result = process_bscan(
                iY, reader, params, Ref_ch1, Ref_ch2, phV, winG, winG_whole,
                nWin, windex, winL, czrg, topLines, SPL, nX, nr, Blength
            )
            
            # 存储结果
            dopu_splitSpectrum[:, :, iY] = result['dopu']
            Strus[:, :, iY] = result['Strus']
            Smap_avg[:, :, :, iY] = result['Smap_avg']
            Smap_rep1[:, :, :, iY] = result['Smap_rep1']
            topLines[:, iY] = result['topLines']
            LA_c_cfg1_avg[:, :, :, iY] = result['LA']
            PhR_c_cfg1_avg[:, :, iY] = result['PhR']
            cumLA_cfg1_avg[:, :, :, iY] = result['cumLA']
            LA_Ms_cfg1_rmBG[:, :, :, iY] = result['LA_raw']
            PhR_Ms_cfg1_rmBG[:, :, iY] = result['PhR_raw']
            cumLA_Ms_cfg1_rmBG[:, :, :, iY] = result['cumLA_raw']
    
    # 保存结果
    if params.tiff.save_dicom:
        print("保存DICOM文件...")
        save_results(
            foutputdir, name, nr, nY, nX, nZcrop, avnum,
            Strus, Smap_avg, Smap_rep1, topLines, dopu_splitSpectrum,
            LA_c_cfg1_avg, PhR_c_cfg1_avg, cumLA_cfg1_avg,
            LA_Ms_cfg1_rmBG, PhR_Ms_cfg1_rmBG, cumLA_Ms_cfg1_rmBG,
            params
        )
    
    # DCM到TIFF转换
    if params.tiff.save_dicom and params.tiff.make_tiff:
        print("\n开始DCM到TIFF转换...")
        try:
            convert_dcm_to_tiff(str(foutputdir), params.tiff.tiff_frame, name)
        except Exception as e:
            print(f"DCM到TIFF转换出错: {e}")
    
    # 显示处理时间
    proc_time = time.time() - file_start_time
    print(f"\n文件 {input_path.name} 处理完成!")
    print(f"输出目录: {foutputdir}")
    print(f"处理时间: {proc_time:.2f} 秒 ({proc_time/60:.2f} 分钟)")


def save_results(
    foutputdir: Path,
    name: str,
    nr: int,
    nY: int,
    nX: int,
    nZcrop: int,
    avnum: int,
    Strus: np.ndarray,
    Smap_avg: np.ndarray,
    Smap_rep1: np.ndarray,
    topLines: np.ndarray,
    dopu_splitSpectrum: np.ndarray,
    LA_c_cfg1_avg: np.ndarray,
    PhR_c_cfg1_avg: np.ndarray,
    cumLA_cfg1_avg: np.ndarray,
    LA_Ms_cfg1_rmBG: np.ndarray,
    PhR_Ms_cfg1_rmBG: np.ndarray,
    cumLA_Ms_cfg1_rmBG: np.ndarray,
    params: ConfigParams
):
    """保存处理结果"""
    from scipy.io import savemat
    
    foutputdir = Path(foutputdir)
    
    # 归一化结构图像
    slice_index = min(100, Strus.shape[2] // 2)
    if slice_index < 1:
        slice_index = 0
    SS = Strus[:, :, slice_index]
    strUrg = SS.max() - 5
    strLrg = SS.min() + 5
    Struc = (Strus - strLrg) / (strUrg - strLrg + 1e-10)
    Struc = np.clip(Struc, 0, 1)
    
    # 保存原始结构图像
    dicom_write((255 * Struc).astype(np.uint8), str(foutputdir / f"{name}_1-1_Struc.dcm"))
    
    # 创建带边界的结构图像
    Struc_with_boundary = Struc.copy()
    for iY in range(Struc.shape[2]):
        for iX in range(Struc.shape[1]):
            surface_pos = int(round(topLines[iX, iY]))
            if 1 < surface_pos <= Struc.shape[0]:
                Struc_with_boundary[:surface_pos-1, iX, iY] = 0
    
    dicom_write((255 * Struc_with_boundary).astype(np.uint8), str(foutputdir / f"{name}_1-2_Struc_with_boundary.dcm"))
    dicom_write((255 * (Smap_rep1 / 2 + 0.5)).astype(np.uint8), str(foutputdir / f"{name}_1-3_1rep-Stokes.dcm"))
    dicom_write((255 * (Smap_avg / 2 + 0.5)).astype(np.uint8), str(foutputdir / f"{name}_1-3_4rep-Stokes.dcm"))
    
    # DOPU处理
    dopu_thresholded = dopu_splitSpectrum.copy()
    dopu_thresholded[dopu_thresholded <= 0.55] = 0
    
    dopu_with_boundary = dopu_splitSpectrum.copy()
    for iY in range(dopu_with_boundary.shape[2]):
        for iX in range(dopu_with_boundary.shape[1]):
            surface_pos = int(round(topLines[iX, iY]))
            if 1 < surface_pos <= dopu_with_boundary.shape[0]:
                dopu_with_boundary[:surface_pos-1, iX, iY] = 0
    
    dicom_write((255 * dopu_thresholded).astype(np.uint8), str(foutputdir / f"{name}_1-4_dopu_SS.dcm"))
    dicom_write((255 * dopu_with_boundary).astype(np.uint8), str(foutputdir / f"{name}_1-5_dopu_SS_with_boundary.dcm"))
    
    # 保存分割结果
    if not params.processing.has_seg:
        czrg = np.arange(1, nZcrop + 1)
        savemat(str(foutputdir / f"{name}.mat"), {'topLines': topLines, 'czrg': czrg})
    
    # 保存偏振参数
    rotAngle = 440
    PRRrg = [0, 0.5]
    np.savetxt(str(foutputdir / f"{name}_2-0_PhRRg.txt"), PRRrg)
    
    parula = parula_colormap(256)
    
    # 处理每个B-scan的颜色编码
    PRRc = np.zeros((nZcrop - avnum, nX, 3, nY), dtype=np.uint8)
    cumLA_cfg_hsv = np.zeros_like(PRRc, dtype=float)
    LA_cfg_hsv = np.zeros_like(PRRc, dtype=float)
    PRRc_rmBG = np.zeros_like(PRRc)
    cumLA_Ms_cfg1_rmBG_hsv = np.zeros_like(PRRc, dtype=float)
    LA_Ms_cfg1_rmBG_hsv = np.zeros_like(PRRc, dtype=float)
    
    for iY in range(nY):
        # 相位延迟颜色编码
        phr_gray = mat2gray(PhR_c_cfg1_avg[:, :, iY], PRRrg)
        phr_idx = (phr_gray * 255).astype(np.uint8)
        PRRc[:, :, :, iY] = (parula[phr_idx] * 255).astype(np.uint8)
        
        # 累积双折射HSV颜色编码
        cumLA_cfg_hsv[:, :, :, iY] = qu_coloring(cumLA_cfg1_avg[:, :, :, iY], rotAngle)
        LA_cfg_hsv[:, :, :, iY] = qu_coloring(LA_c_cfg1_avg[:, :, :, iY], rotAngle)
        
        # 去背景版本
        phr_gray_rmBG = mat2gray(PhR_Ms_cfg1_rmBG[:, :, iY], PRRrg)
        phr_idx_rmBG = (phr_gray_rmBG * 255).astype(np.uint8)
        PRRc_rmBG[:, :, :, iY] = (parula[phr_idx_rmBG] * 255).astype(np.uint8)
        
        cumLA_Ms_cfg1_rmBG_hsv[:, :, :, iY] = qu_coloring(cumLA_Ms_cfg1_rmBG[:, :, :, iY], rotAngle)
        LA_Ms_cfg1_rmBG_hsv[:, :, :, iY] = qu_coloring(LA_Ms_cfg1_rmBG[:, :, :, iY], rotAngle)
    
    # 保存偏振参数DICOM文件
    dicom_write((255 * (cumLA_cfg1_avg / 2 + 0.5)).astype(np.uint8), str(foutputdir / f"{name}_2-1_cumLA-cfg1-{nr}repAvg.dcm"))
    dicom_write((255 * cumLA_cfg_hsv).astype(np.uint8), str(foutputdir / f"{name}_2-2_cumLA-cfg1-{nr}repAvg_hsvColoring.dcm"))
    dicom_write((255 * (LA_c_cfg1_avg / 2 + 0.5)).astype(np.uint8), str(foutputdir / f"{name}_2-3_drLA-cfg1-{nr}repAvg.dcm"))
    dicom_write((255 * LA_cfg_hsv).astype(np.uint8), str(foutputdir / f"{name}_2-4_drLA-cfg1-{nr}repAvg_hsvColoring.dcm"))
    dicom_write(PRRc, str(foutputdir / f"{name}_2-5_PhR-cfg1-{nr}repAvg.dcm"))
    
    dicom_write((255 * (cumLA_Ms_cfg1_rmBG / 2 + 0.5)).astype(np.uint8), str(foutputdir / f"{name}_2-6_cumLA_rmBG-cfg1-{nr}repAvg.dcm"))
    dicom_write((255 * cumLA_Ms_cfg1_rmBG_hsv).astype(np.uint8), str(foutputdir / f"{name}_2-7_cumLA_rmBG-cfg1-{nr}repAvg_hsvColoring.dcm"))
    dicom_write((255 * (LA_Ms_cfg1_rmBG / 2 + 0.5)).astype(np.uint8), str(foutputdir / f"{name}_2-8_drLA_rmBG-cfg1-{nr}repAvg.dcm"))
    dicom_write((255 * LA_Ms_cfg1_rmBG_hsv).astype(np.uint8), str(foutputdir / f"{name}_2-9_drLA_rmBG-cfg1-{nr}repAvg_hsvColoring.dcm"))
    dicom_write(PRRc_rmBG, str(foutputdir / f"{name}_2-10_PhR_rmBG-cfg1-{nr}repAvg.dcm"))
    
    print(f"已保存结果到: {foutputdir}")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else "./output"
        
        process_single_file(input_file, output_dir)
    else:
        print("用法: python psoct_processor.py <oct文件路径> [输出目录]")
