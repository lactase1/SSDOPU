%% ====================================================================================
% 函数名: calLAPhRALL (完整版偏振参数计算)
% 功能: 计算完整的偏振参数，包括处理前和处理后的结果
% 输入参数:
%   IMG_ch1, IMG_ch2 - 双通道复数OCT信号
%   test_seg_top - 组织表面位置
%   dopu_splitSpec_M - 分光谱DOPU矩阵
%   kRL, kRU - 滤波核范围参数
%   h1, h2 - 多尺度高斯滤波核
%   Avnum - 平均层数
%   wovWinF - 滤波模式(可选参数，默认为0)
%   enableDopuPhaseSupp - DOPU相位抑制开关(可选，默认为1)
%   verbose - 详细输出等级(可选，默认为0: 静默)
%   enableOutputAdaptive - 启用输出端混合滤波(可选，默认为0)
%   kRL_output - 输出滤波核下限(可选，默认为13)
%   kRU_output - 输出滤波核上限(可选，默认为21)
%   outputDopuThreshold - 输出滤泥DOPU阈值(可选，默认为0.4)
%   enableBottomLayerPhaseReduction - 启用底层相位延迟减小(可选，默认为0)
%   bottomLayerDepth - 底层深度范围(可选，默认为80)
% 输出参数:
%   LA - 处理后的局部双折射 [(Z-Avnum)×X×3]
%   PhR - 处理后的相位延迟 [(Z-Avnum)×X]
%   cumLA - 处理后的累积双折射 [(Z-Avnum)×X×3]
%   LA_raw - 原始局部双折射(去背景前)
%   PhR_raw - 原始相位延迟(去背景前)
%   cumLA_raw - 原始累积双折射(去背景前)
% 特点:
%   - 提供完整的处理流程，包括原始和处理后的结果
%   - 支持两种滤波模式的完整实现
%   - 支持输出端自适应滤波，根据DOPU动态调整光轴/延迟滤波
%   - 支持底层相位延迟减小，降低深层区域噪声
%   - 用于深入的偏振分析和算法验证
% ====================================================================================
function [LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = calLAPhRALL(IMG_ch1,IMG_ch2,test_seg_top,dopu_splitSpec_M,kRL,kRU,h1,h2,Avnum,wovWinF,enableDopuPhaseSupp, verbose, enableOutputAdaptive, kRL_output, kRU_output, outputDopuThreshold, enableBottomLayerPhaseReduction, bottomLayerDepth)
    % 处理可选参数
    if nargin < 10, wovWinF = 0; end  % 直接设置wovWinF而不是params.mode.wovWinF
    if nargin < 11, enableDopuPhaseSupp = 1; end
    if nargin < 12, verbose = 0; end  % 0 = 静默, 1 = 简短, 2 = 详细
    if nargin < 13, enableOutputAdaptive = 0; end  % 默认禁用自适应输出滤泥
    if nargin < 14, kRL_output = 13; end  % 默认滤泥核下限
    if nargin < 15, kRU_output = 21; end  % 默认滤泥核上限
    if nargin < 16, outputDopuThreshold = 0.4; end  % 默认输出滤泥DOPU阈值
    if nargin < 17, enableBottomLayerPhaseReduction = 0; end  % 默认禁用底层相位延迟减小
    if nargin < 18, bottomLayerDepth = 80; end  % 默认底层深度80层

    % 获取数据维度
    [nZ, nX] = size(IMG_ch1);

    % 初始化处理后的输出矩阵
    LA = zeros(nZ-Avnum, nX, 3);       % 使用传入的Avnum参数
    PhR = zeros(nZ-Avnum, nX);         % 使用传入的Avnum参数
    cumLA = LA;                        % 处理后累积双折射

    % 初始化原始输出矩阵
    LA_raw = zeros(nZ-Avnum, nX, 3);   % 使用传入的Avnum参数
    PhR_raw = zeros(nZ-Avnum, nX);     % 使用传入的Avnum参数
    cumLA_raw = LA;                    % 原始累积双折射

    % 计算Stokes参数
    [ES0, ES1, ES2, ES3] = cumulativeQUV(IMG_ch1, IMG_ch2);

    % 信号强度检查
    if sum(ES0(:)) < 5, return; end

    % 计算归一化Stokes分量
    EQm = ES1 ./ ES0;   % Q分量(水平-垂直偏振差)
    EUm = ES2 ./ ES0;   % U分量(±45°偏振差)
    EVm = ES3 ./ ES0;   % V分量(左右旋圆偏振差)

    % 初始化结构矩阵
    Stru_E = zeros(nZ, nX);

    % 根据滤波模式进行预处理
    if wovWinF == 1
        % 模式1: 固定高斯滤波
        EQmm = imfilter(EQm, h1, 'replicate');
        EUmm = imfilter(EUm, h1, 'replicate');
        EVmm = imfilter(EVm, h1, 'replicate');
    else
        % 模式0: 自适应DOPU滤波
        [EQmm] = vWinAvgFiltOpt(EQm, dopu_splitSpec_M, kRL, kRU);
        [EUmm] = vWinAvgFiltOpt(EUm, dopu_splitSpec_M, kRL, kRU);
        [EVmm] = vWinAvgFiltOpt(EVm, dopu_splitSpec_M, kRL, kRU);
    end

    % 调用核心算法，同时获得处理前后的完整结果
    [LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw] = FreeSpace_PSOCT_Optimized(EQmm, EUmm, EVmm, Stru_E, test_seg_top, h1, h2, Avnum, dopu_splitSpec_M, enableDopuPhaseSupp, 0.42, enableOutputAdaptive, kRL_output, kRU_output, outputDopuThreshold, enableBottomLayerPhaseReduction, bottomLayerDepth);
end
