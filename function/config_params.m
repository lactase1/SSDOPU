function params = config_params()
%% ====================================================================================
% 函数名: config_params
% 功能: PS-OCT处理参数配置文件
% 输出: params - 包含所有处理参数的结构体
% 
% 使用方法:
%   params = config_params();
%   make_tiff = params.tiff.make_tiff;
%   
% 更新时间: 2025年9月17日
% 说明: 集中管理所有PS-OCT处理参数，便于维护和调整
% ====================================================================================

%% 初始化参数结构体
params = struct();

%% ========== TIFF生成控制参数 ==========
params.tiff.make_tiff = 1;        % 1: 生成TIFF文件; 0: 不生成
params.tiff.tiff_frame = 270;     % 要提取的帧号(默认160，即中间帧)
params.tiff.saveDicom = 1;        % 是否保存DICOM文件 (1:保存, 0:不保存)

%% ========== DOPU增强处理参数 ==========
params.dopu.enable_enhanced = 1;               % 是否启用DOPU增强处理 (1:启用, 0:禁用)
params.dopu.noise_threshold = 0.2;             % DOPU噪声阈值(提高以更好去除噪声，改善分层清晰度)
params.dopu.contrast_enhance = 0.6;            % DOPU对比度增强因子(降低以进一步增强对比度)
params.dopu.create_overlay = 1;                % 是否创建结构图像叠加DOPU (1:启用, 0:禁用)

%% ========== 基础处理参数 ==========
params.processing.disp_coef = -20.1;           % 色散补偿系数
params.processing.do_PhComp = 1;               % 是否进行相位补偿 (1:开启, 0:关闭)
params.processing.do_medianshift = 1;          % 是否进行中值偏移校正
params.processing.do_reg = 0;                  % 是否进行帧间配准(运动去除)
params.processing.useref = 1;                  % 参考信号模式 (1:使用前50k A-line, 0:使用背景, -1:不使用)
params.processing.show_img = 0;                % 是否显示中间结果图像
params.processing.iy = 1;                      % Y方向步长(通常为1)
params.processing.hasSeg = 0;                  % 是否已有分割结果(.mat文件)
params.range.setZrg = 0;

%% ========== 并行处理设置 ==========
params.parallel.LocalUseMpiexec = false;       % 并行处理MPI设置

%% ========== 处理模式配置 ==========
% 配置1: 标准双折射分析 (LA calculation)
params.mode.do_cfg1 = 1;                       % 是否启用配置1
% 配置2: 相位延迟分析 (Phase Retardation calculation) 
params.mode.do_cfg2 = 0;                       % 是否启用配置2
% 滤波模式
params.mode.wovWinF = 1;                       % 滤波模式 (1:固定高斯滤波, 0:自适应DOPU滤波)

%% ========== 分光谱DOPU设置 ==========
params.dopu.do_ssdopu = 1;                     % 是否启用分光谱DOPU (1:启用, 0:禁用)
params.dopu.do_avg = 1;                        % 是否使用平均方法 (1:启用, 0:禁用)
params.dopu.do_eig = 0;                        % 是否使用特征值方法 (1:启用, 0:禁用)

%% ========== 偏振分析参数 ==========
% 【优化说明】基于PS-OCT技术特性调整偏振分析参数：
% - 巩膜组织具有显著双折射特性，需要精确的相位测量
% - Avnum影响信号平滑度与空间分辨率的平衡
% - 滤波核范围影响DOPU计算的稳定性和精度

% 平均层数设置
params.polarization.Avnum_initial = 7;         % 初始平均层数(适当减小以提高分辨率)
params.polarization.Avnum = 4;                 % DDG测试用平均层数(统一使用以保持一致性)

% 配置1滤波核范围 (用于局部双折射LA计算)
% 调整为更适合巩膜结构特征的范围
params.polarization.kRL_cfg1 = 8;              % 配置1滤波核下限(减小以保持细节)
params.polarization.kRU_cfg1 = 15;              % 配置1滤波核上限(适当增加覆盖范围)

% 配置2滤波核范围 (用于相位延迟PhR计算)
% 相位延迟测量需要更大的统计区域以提高精度  
params.polarization.kRL_cfg2 = 2;              % 配置2滤波核下限(微调下限)
params.polarization.kRU_cfg2 = 15;             % 配置2滤波核上限(减小以平衡精度与分辨率)

%% ========== 高斯滤波核设置 ==========
% 【优化说明】基于PS-OCT巩膜检测技术需求的参数调整：
% - 巩膜检测需要保持细节结构的同时抑制散斑噪声
% - 双折射测量要求在噪声抑制与空间分辨率间平衡
% - 考虑人眼巩膜组织的结构特征尺度
% - 基于DOPU图像质量反馈，增大滤波核以改善分层清晰度

% 小尺度高斯核 (用于细节保持和初步降噪) - 优化版
params.filters.h1_size = [7 7];                % 高斯核1尺寸 (增大以更好抑制噪声)
params.filters.h1_sigma = 4;                 % 高斯核1标准差 (提高以改善分层效果)
params.filters.h1 = fspecial('gaussian', params.filters.h1_size, params.filters.h1_sigma);

% 中尺度高斯核 (用于结构增强和背景平滑) - 优化版
params.filters.h2_size = [15 15];              % 高斯核2尺寸 (适当增大以增强结构连续性)
params.filters.h2_sigma = 6.0;                 % 高斯核2标准差 (提高以更好平滑背景)
params.filters.h2 = fspecial('gaussian', params.filters.h2_size, params.filters.h2_sigma);

%% ========== 数据处理范围设置 ==========
% 【优化说明】根据眼部结构和巩膜位置特点调整处理范围
% 注意: 实际范围将根据数据文件的SPL自动调整，确保不超过SPL/2
params.range.czrg_start = 15;                  % Z方向处理起始位置(略提前以包含更多巩膜信息)
params.range.czrg_end = 280;                   % Z方向处理结束位置(保守设置，避免超出数据范围)
params.range.setZrg = 0;                       % 深度范围限制(0:处理全部深度以获得完整巩膜信息)
params.range.Thr = 150;                        % OCT信号阈值(略降低以包含更多弱信号区域)

%% ========== 分光谱DOPU窗口设置 ==========
% 【优化说明】调整分光谱参数以优化巩膜检测性能
params.splitspec.nWin = 11;                    % 分光谱窗口数量(增加以提高频谱分辨率)
params.splitspec.tukey_alpha = 0.15;           % Tukey窗口参数(减小以减少频谱泄漏)

%% ========== 表面检测参数 ==========
% 【优化说明】针对眼部结构特点调整表面检测参数  
params.surface.surf_threshold = 0.20;          % 表面检测阈值(略降低以检测更多细微结构)
params.surface.median_window = 7;              % 中值滤波窗口大小(减小以保持边界清晰)
params.surface.smooth_window = 9;              % 平滑滤波窗口大小(减小以保持局部特征)

%% ========== 验证参数合理性 ==========
% 【增强验证】检查所有关键参数的合理性
% 检查滤波核范围参数
if params.polarization.kRL_cfg1 >= params.polarization.kRU_cfg1
    warning('配置1滤波核范围设置有误: kRL_cfg1应小于kRU_cfg1');
end

if params.polarization.kRL_cfg2 >= params.polarization.kRU_cfg2
    warning('配置2滤波核范围设置有误: kRL_cfg2应小于kRU_cfg2');
end

% 检查处理范围参数
if params.range.czrg_start >= params.range.czrg_end
    warning('Z方向处理范围设置有误: 起始位置应小于结束位置');
end

% 检查滤波核尺寸合理性
if any(params.filters.h1_size < 1) || any(params.filters.h2_size < 1)
    warning('滤波核尺寸设置有误: 尺寸应为正整数');
end

% 检查sigma参数合理性
if params.filters.h1_sigma <= 0 || params.filters.h2_sigma <= 0
    warning('高斯核标准差设置有误: 标准差应为正数');
end

%% ========== 参数信息输出(可选) ==========
if nargout == 0
    % 如果直接调用函数不接收输出，则显示参数摘要
    fprintf('\n=== PS-OCT处理参数配置摘要 ===\n');
    fprintf('TIFF生成: %s (帧数: %d)\n', ...
        iif(params.tiff.make_tiff, '启用', '禁用'), params.tiff.tiff_frame);
    fprintf('DICOM保存: %s\n', iif(params.tiff.saveDicom, '启用', '禁用'));
    fprintf('分光谱DOPU: %s\n', iif(params.dopu.do_ssdopu, '启用', '禁用'));
    fprintf('DOPU增强处理: %s (噪声阈值: %.2f, 对比度因子: %.2f)\n', ...
        iif(params.dopu.enable_enhanced, '启用', '禁用'), params.dopu.noise_threshold, params.dopu.contrast_enhance);
    fprintf('DOPU结构叠加: %s\n', iif(params.dopu.create_overlay, '启用', '禁用'));
    fprintf('配置1(LA): %s, 配置2(PhR): %s\n', ...
        iif(params.mode.do_cfg1, '启用', '禁用'), iif(params.mode.do_cfg2, '启用', '禁用'));
    fprintf('滤波模式: %s\n', iif(params.mode.wovWinF, '固定高斯', '自适应DOPU'));
    fprintf('平均层数: %d\n', params.polarization.Avnum);
    fprintf('==========================\n\n');
end