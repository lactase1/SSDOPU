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
params.tiff.tiff_frame = 35;     % 要提取的帧号(默认160，即中间帧)
params.tiff.saveDicom = 1;        % 是否保存DICOM文件 (1:保存, 0:不保存)

%% ========== 基础处理参数 ==========
params.processing.disp_coef = -20.1;           % 色散补偿系数
params.processing.do_PhComp = 1;               % 是否进行相位补偿 (1:开启, 0:关闭)
params.processing.do_medianshift = 1;          % 是否进行中值偏移校正
params.processing.do_reg = 0;                  % 是否进行帧间配准(运动去除)
params.processing.useref = 1;                  % 参考信号模式 (1:使用前50k A-line, 0:使用背景, -1:不使用)
params.processing.show_img = 0;                % 是否显示中间结果图像
params.processing.iy = 1;                      % Y方向步长(通常为1)
params.processing.hasSeg = 1;                  % 是否已有分割结果(.mat文件)
params.processing.max_frames = 35;              % 最大处理帧数 (0:处理所有帧, >0:限制帧数)
params.range.setZrg = 0;

%% ========== 并行处理设置 ==========
params.parallel.LocalUseMpiexec = false;       % 并行处理MPI设置
% 可选：限制并行池最大 worker 数 (0 表示由代码自动选择)
params.parallel.maxWorkers = 48;

%% ========== 处理模式配置 ==========
% 处理模式配置（当前仅保留滤波模式选项）
params.mode.wovWinF = 0;                       % 滤波模式 (1:固定高斯滤波, 0:自适应DOPU滤波)

%% ========== 分光谱DOPU设置 ==========
params.dopu.do_ssdopu = 1;                     % 是否启用分光谱DOPU (1:启用, 0:禁用)
% 分裂谱模式: 'overlap9' (默认，9 段半重叠) 或 'nonoverlap5' (5 段无重叠)
params.dopu.ss_mode = 'overlap9';
% 窗口 Tukey alpha（可由 nonoverlap5 函数使用）
params.dopu.win_alpha = 0.25;

%% ========== 空间DOPU设置 ==========
params.dopu.do_spatial = 0;                    % 是否启用空间DOPU (1:启用, 0:禁用)
% 空间DOPU基于3x3邻域平均计算

%% ========== 组合DOPU设置 ==========
params.dopu.do_combined = 1;                   % 是否启用组合DOPU (分裂谱+空间) (1:启用, 0:禁用)
% 组合DOPU先计算分裂谱DOPU，再进行空间滤波

%% ========== 偏振分析参数 ==========
% 【优化说明】基于PS-OCT技术特性调整偏振分析参数：
% - 巩膜组织具有显著双折射特性，需要精确的相位测量
% - Avnum影响信号平滑度与空间分辨率的平衡
% - 滤波核范围影响DOPU计算的稳定性和精度

% 平均层数设置
params.polarization.Avnum = 3;                 % DDG测试用平均层数(统一使用以保持一致性)
params.polarization.enableDopuPhaseSupp = 0;  % 1: 使用DOPU自适应相位抑制; 0: 关闭该功能

% 配置1滤波核范围 (用于局部双折射LA计算)
% 调整为更适合巩膜结构特征的范围
params.polarization.kRL_cfg1 = 3;              % 配置1滤波核下限(减小以保持细节)
params.polarization.kRU_cfg1 = 21;              % 配置1滤波核上限(适当增加覆盖范围)

% （已移除配置2相关参数，当前仅使用配置1的滤波范围）

%% ========== 高斯滤波核设置 ==========
% 【优化说明】基于PS-OCT巩膜检测技术需求的参数调整：
% - 巩膜检测需要保持细节结构的同时抑制散斑噪声
% - 双折射测量要求在噪声抑制与空间分辨率间平衡
% - 考虑人眼巩膜组织的结构特征尺度
% - 基于DOPU图像质量反馈，增大滤波核以改善分层清晰度

% 小尺度高斯核 (用于细节保持和初步降噪) - 细节优先
params.filters.h1_size = [3 3];                % 高斯核1尺寸 (减小以增强局部纹理)
params.filters.h1_sigma = 1.2;                % 高斯核1标准差 (降低以避免过度平滑)
params.filters.h1 = fspecial('gaussian', params.filters.h1_size, params.filters.h1_sigma);

% 中尺度高斯核 (用于结构增强和背景平滑) - 背景优先
params.filters.h2_size = [20 20];               % 高斯核2尺寸 (显著增大以平滑大尺度背景)
params.filters.h2_sigma = 4;                  % 高斯核2标准差 (增大以抑制深层噪声)
params.filters.h2 = fspecial('gaussian', params.filters.h2_size, params.filters.h2_sigma);

% （已删除若干未使用的注释/备用参数，以保持配置简洁）

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

% （配置2 已弃用，相关检查移除）

% % 检查处理范围参数
% if params.range.czrg_start >= params.range.czrg_end
%     warning('Z方向处理范围设置有误: 起始位置应小于结束位置');
% end

% 检查滤波核尺寸合理性
if any(params.filters.h1_size < 1) || any(params.filters.h2_size < 1)
    warning('滤波核尺寸设置有误: 尺寸应为正整数');
end

% % 检查sigma参数合理性
% if params.filters.h1_sigma <= 0 || params.filters.h2_sigma <= 0
%     warning('高斯核标准差设置有误: 标准差应为正数');
% end

%% ========== 参数信息输出(可选) ==========
if nargout == 0
    % 如果直接调用函数不接收输出，则显示参数摘要
    fprintf('\n=== PS-OCT处理参数配置摘要 ===\n');
    fprintf('TIFF生成: %s (帧: %d)\n', iif(params.tiff.make_tiff, '启用', '禁用'), params.tiff.tiff_frame);
    fprintf('DICOM保存: %s\n', iif(params.tiff.saveDicom, '启用', '禁用'));
    fprintf('分光谱DOPU: %s (模式: %s)\n', iif(params.dopu.do_ssdopu, '启用', '禁用'), params.dopu.ss_mode);
    fprintf('空间DOPU: %s\n', iif(params.dopu.do_spatial, '启用', '禁用'));
    fprintf('组合DOPU: %s\n', iif(params.dopu.do_combined, '启用', '禁用'));
    fprintf('DOPU计算方法: avg (固定)\n');
    fprintf('滤波模式: %s\n', iif(params.mode.wovWinF, '固定高斯', '自适应DOPU'));
    fprintf('平均层数 (Avnum): %d\n', params.polarization.Avnum);
        fprintf('DOPU相位抑制: %s\n', iif(params.polarization.enableDopuPhaseSupp, '启用', '禁用'));
    fprintf('==========================\n\n');
end