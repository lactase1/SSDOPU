%% ====================================================================================
% 分裂谱DOPU计算函数使用指南
% 函数名: calculateSplitSpectrumDOPU.m
% 版本: 1.0
% 更新日期: 2024年10月
% ====================================================================================

%% 1. 函数概述
% calculateSplitSpectrumDOPU 函数用于计算分裂谱偏振均匀度(Degree of Polarization Uniformity)
% 该函数将OCT信号的频谱分割成多个子频带，对每个子频带分别计算偏振参数，
% 通过频谱平均提高信噪比和空间分辨率。

%% 2. 算法原理
% - 将OCT信号频谱分割成多个重叠的子频带 (默认9个窗口)
% - 对每个子频带应用高斯窗口并进行FFT变换
% - 计算每个子频带的Stokes参数 (S0, S1, S2, S3)
% - 对所有子频带进行平均
% - 计算DOPU = sqrt(<S1/S0>² + <S2/S0>² + <S3/S0>²)

%% 3. 函数签名
% [dopu_splitSpectrum, dopu_ss] = calculateSplitSpectrumDOPU(...
%     Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg)

%% 4. 输入参数详细说明

%% 4.1 信号数据
% Bd1 - 通道1的OCT信号矩阵
%       类型: double
%       尺寸: [信号长度 × X像素 × 重复次数]
%       说明: 经过背景减除和色散校正的OCT信号
%
% Bd2 - 通道2的OCT信号矩阵
%       类型: double
%       尺寸: [信号长度 × X像素 × 重复次数]
%       说明: 第二个偏振通道的OCT信号

%% 4.2 参数结构体
% params - 参数结构体，必须包含以下字段:
%   params.dopu.do_ssdopu - 是否进行分裂谱计算
%       类型: logical
%       默认值: true (进行计算) 或 false (返回默认值1)

%% 4.3 信号尺寸参数
% SPL - 信号长度 (FFT点数)
%       类型: integer
%       典型值: 1024, 2048等
%
% nX - X方向像素数 (A-scan数量)
%       类型: integer
%       典型值: 512, 1024等
%
% nr - 重复次数 (repetitions per B-scan)
%       类型: integer
%       典型值: 1, 2, 4等

%% 4.4 分裂谱参数
% nWin - 分裂谱窗口数
%       类型: integer
%       默认值: 9
%       说明: 频谱被分割的子频带数量
%
% windex - 窗口起始索引数组
%       类型: integer array
%       尺寸: [1 × nWin]
%       说明: 每个窗口在信号中的起始位置
%       示例: windex = 1 : winL/2 : SPL;
%
% winL - 单个窗口长度
%       类型: integer
%       计算公式: winL = 2*SPL/(nWin+1)
%       说明: 每个子频带的长度
%
% winG - 高斯窗口函数
%       类型: double array
%       尺寸: [winL × 1]
%       说明: Tukey窗口 (余弦锥形窗口)
%       生成方法: winG = tukeywin(winL, 0.25);

%% 4.5 Z方向裁剪参数
% czrg - Z方向裁剪范围
%       类型: integer array
%       示例: czrg = 1:320 或 czrg = 60:380
%       说明: 选择感兴趣的深度范围

%% 5. 输出参数详细说明

%% 5.1 dopu_splitSpectrum
% 类型: double
% 尺寸: [nZcrop × nX] 其中 nZcrop = length(czrg)
% 说明: 最终的分裂谱DOPU结果，对所有重复次数进行了平均
% 取值范围: [0, 1]，1表示完全偏振，0表示非偏振

%% 5.2 dopu_ss
% 类型: double
% 尺寸: [nZcrop × nX × nr]
% 说明: 各重复次数的分裂谱DOPU结果
% 用途: 用于分析重复测量的一致性

%% 6. 调用示例

%% 6.1 基本调用 (在PS-OCT处理流程中)
% 假设已经在主函数中定义了所有必要参数
[dopu_result, dopu_all_reps] = calculateSplitSpectrumDOPU(...
    Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);

%% 6.2 简化调用 (不进行分裂谱计算)
params.dopu.do_ssdopu = false;
[dopu_result, ~] = calculateSplitSpectrumDOPU(...
    Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);
% 返回 dopu_result = ones(nZcrop, nX)

%% 6.3 完整参数准备示例
% % 基本参数设置
% SPL = 1024;           % FFT点数
% nX = 512;             % A-scan数量
% nr = 4;               % 重复次数
% nWin = 9;             % 分裂谱窗口数
%
% % 计算窗口参数
% winL = 2*SPL/(nWin+1);                    % 窗口长度
% winG = tukeywin(winL, 0.25);              % 高斯窗口
% windex = 1 : winL/2 : SPL;                % 窗口起始位置
%
% % Z方向范围
% czrg = 1:320;         % 感兴趣的深度范围
%
% % 调用函数
% [dopu_splitSpectrum, dopu_ss] = calculateSplitSpectrumDOPU(...
%     Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);

%% 7. 错误处理
% 函数包含以下错误检查:
% - 参数数量验证
% - 输入信号非空检查
% - 信号尺寸一致性检查
% - 尺寸与参数匹配检查

%% 8. 性能考虑
% - 计算复杂度: O(nZcrop × nX × nr × nWin × log(SPL))
% - 内存使用: 主要消耗在FFT结果存储 Bimg1, Bimg2
% - 建议: 对于大矩阵，考虑分块处理

%% 9. 依赖关系
% 依赖函数:
% - cumulativeQUV: 计算Stokes参数
% MATLAB工具箱:
% - Signal Processing Toolbox (tukeywin, fft)

%% 10. 版本历史
% v1.0 (2024-10) - 初始版本，从主处理脚本中抽离
% - 支持标准分裂谱DOPU计算
% - 包含详细的参数验证
% - 提供完整的错误处理

%% 11. 参考文献
% [1] 分裂谱OCT技术用于提高空间分辨率
% [2] 偏振敏感OCT中的DOPU计算方法
% ====================================================================================