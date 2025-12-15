"""
PS-OCT处理参数配置文件
功能: 集中管理所有PS-OCT处理参数
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class TiffParams:
    """TIFF生成控制参数"""
    make_tiff: bool = True          # 是否生成TIFF文件
    tiff_frame: int = 35            # 要提取的帧号
    save_dicom: bool = True         # 是否保存DICOM文件


@dataclass
class ProcessingParams:
    """基础处理参数"""
    disp_coef: float = -20.1        # 色散补偿系数
    do_ph_comp: bool = True         # 是否进行相位补偿
    do_medianshift: bool = True     # 是否进行中值偏移校正
    do_reg: bool = False            # 是否进行帧间配准
    useref: int = 1                 # 参考信号模式 (1:使用前50k A-line, 0:使用背景, -1:不使用)
    show_img: bool = False          # 是否显示中间结果图像
    iy: int = 1                     # Y方向步长
    has_seg: bool = True            # 是否已有分割结果
    max_frames: int = 0             # 最大处理帧数 (0:处理所有帧)


@dataclass
class RangeParams:
    """处理范围参数"""
    set_z_rg: int = 0               # Z方向处理范围 (0:处理全部)


@dataclass
class ParallelParams:
    """并行处理参数"""
    local_use_mpiexec: bool = False
    max_workers: int = 48           # 最大worker数


@dataclass
class ModeParams:
    """处理模式配置"""
    wov_win_f: int = 0              # 滤波模式 (1:固定高斯滤波, 0:自适应DOPU滤波)


@dataclass
class DopuParams:
    """DOPU设置"""
    do_ssdopu: bool = True          # 是否启用分光谱DOPU
    ss_mode: str = 'overlap9'       # 分裂谱模式
    win_alpha: float = 0.25         # 窗口Tukey alpha
    do_spatial: bool = False        # 是否启用空间DOPU
    do_combined: bool = True        # 是否启用组合DOPU
    verbose: bool = False           # 是否打印详细信息


@dataclass
class PolarizationParams:
    """偏振分析参数"""
    avnum: int = 3                  # DDG测试用平均层数
    enable_dopu_phase_supp: bool = False  # 是否使用DOPU自适应相位抑制
    krl_cfg1: int = 3               # 配置1滤波核下限
    kru_cfg1: int = 21              # 配置1滤波核上限


@dataclass
class FilterParams:
    """滤波器参数"""
    h1_size: tuple = (3, 3)         # 高斯核1尺寸
    h1_sigma: float = 1.2           # 高斯核1标准差
    h2_size: tuple = (20, 20)       # 高斯核2尺寸
    h2_sigma: float = 4.0           # 高斯核2标准差
    
    def __post_init__(self):
        """初始化高斯核"""
        self.h1 = self._gaussian_kernel(self.h1_size, self.h1_sigma)
        self.h2 = self._gaussian_kernel(self.h2_size, self.h2_sigma)
    
    @staticmethod
    def _gaussian_kernel(size: tuple, sigma: float) -> np.ndarray:
        """生成高斯核"""
        x = np.arange(size[0]) - (size[0] - 1) / 2
        y = np.arange(size[1]) - (size[1] - 1) / 2
        xx, yy = np.meshgrid(x, y)
        kernel = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
        return kernel / kernel.sum()


@dataclass
class SurfaceParams:
    """表面检测参数"""
    surf_threshold: float = 0.20    # 表面检测阈值
    median_window: int = 7          # 中值滤波窗口大小
    smooth_window: int = 9          # 平滑滤波窗口大小


@dataclass
class ConfigParams:
    """PS-OCT处理参数配置总类"""
    tiff: TiffParams = field(default_factory=TiffParams)
    processing: ProcessingParams = field(default_factory=ProcessingParams)
    range: RangeParams = field(default_factory=RangeParams)
    parallel: ParallelParams = field(default_factory=ParallelParams)
    mode: ModeParams = field(default_factory=ModeParams)
    dopu: DopuParams = field(default_factory=DopuParams)
    polarization: PolarizationParams = field(default_factory=PolarizationParams)
    filters: FilterParams = field(default_factory=FilterParams)
    surface: SurfaceParams = field(default_factory=SurfaceParams)
    
    def validate(self):
        """验证参数合理性"""
        if self.polarization.krl_cfg1 >= self.polarization.kru_cfg1:
            print("警告: 配置1滤波核范围设置有误: krl_cfg1应小于kru_cfg1")
        
        if any(s < 1 for s in self.filters.h1_size) or any(s < 1 for s in self.filters.h2_size):
            print("警告: 滤波核尺寸设置有误: 尺寸应为正整数")
    
    def print_summary(self):
        """打印参数摘要"""
        print("\n=== PS-OCT处理参数配置摘要 ===")
        print(f"TIFF生成: {'启用' if self.tiff.make_tiff else '禁用'} (帧: {self.tiff.tiff_frame})")
        print(f"DICOM保存: {'启用' if self.tiff.save_dicom else '禁用'}")
        print(f"分光谱DOPU: {'启用' if self.dopu.do_ssdopu else '禁用'} (模式: {self.dopu.ss_mode})")
        print(f"空间DOPU: {'启用' if self.dopu.do_spatial else '禁用'}")
        print(f"组合DOPU: {'启用' if self.dopu.do_combined else '禁用'}")
        print(f"滤波模式: {'固定高斯' if self.mode.wov_win_f else '自适应DOPU'}")
        print(f"平均层数 (Avnum): {self.polarization.avnum}")
        print(f"DOPU相位抑制: {'启用' if self.polarization.enable_dopu_phase_supp else '禁用'}")
        print("==========================\n")


def config_params() -> ConfigParams:
    """
    创建并返回默认配置参数
    
    Returns:
        ConfigParams: 包含所有处理参数的配置对象
    """
    params = ConfigParams()
    params.validate()
    return params


if __name__ == "__main__":
    # 测试配置
    params = config_params()
    params.print_summary()
