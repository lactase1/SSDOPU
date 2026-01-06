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
params.tiff.make_tiff = 0;        % 1: 生成TIFF文件; 0: 不生成
params.tiff.tiff_frame = 35;     % 要提取的帧号(默认160，即中间帧)
params.tiff.saveDicom = 1;        % 是否保存DICOM文件 (1:保存, 0:不保存)

%% ========== 基础处理参数 ==========
params.processing.disp_coef = -20.1;           % 色散补偿系数
params.processing.do_PhComp = 1;               % 是否进行相位/色散补偿 (1:启用, 0:禁用)
params.processing.do_medianshift = 1;          % 是否进行中值偏移校正
params.processing.do_reg = 0;                  % 是否进行帧间配准(运动去除)
params.processing.useref = 1;                  % 参考信号模式 (1:使用前50k A-line, 0:使用背景, -1:不使用)
params.processing.show_img = 0;                % 是否显示中间结果图像
params.processing.iy = 1;                      % Y方向步长(通常为1)
params.processing.hasSeg = 1;                  % 是否已有分割结果(.mat文件)
params.processing.enable_flatten_enface = 0;   % 1: 启用展平并保存展平体 & 生成 En-face, 0: 禁用
params.processing.enable_enface_noflat = 0;    % 1: 生成非展平En-face切片（直接从原始数据切片）, 0: 禁用
params.processing.max_frames = 0;              % 最大处理帧数 (0:处理所有帧, >0:限制帧数)
params.range.setZrg = 0;
params.parallel.batchSize = 500;

%% ========== 并行处理设置 ==========
params.parallel.LocalUseMpiexec = false;       % 并行处理MPI设置
% 可选：限制并行池最大 worker 数 (0 表示由代码自动选择)
params.parallel.maxWorkers = 47;
% 并行相关资源限制: 最大可用内存 (GB) 用于一次性预加载阈值
params.parallel.maxMemGB = 100;
% 是否在函数结束时自动关闭并行池 (false = 保持池存活以节省启动时间)
params.parallel.autoClosePool = false;

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

% 平均层数设置(MAX_AVNUM = 19)
params.polarization.Avnum = 3;                 % DDG测试用平均层数(统一使用以保持一致性)
params.polarization.enableDopuPhaseSupp = 0;  % 1: 使用DOPU自适应相位抑制; 0: 关闭该功能

% 配置1滤波核范围 (用于局部双折射LA计算)
% 调整为更适合巩膜结构特征的范围
params.polarization.kRL_cfg1 = 3;              % 配置1滤波核下限(减小以保持细节)
params.polarization.kRU_cfg1 = 21;              % 配置1滤波核上限(适当增加覆盖范围)

% 相位延迟范围 (用于PhR彩色编码)
params.polarization.PRRrg = [0.06 0.35];          % PhR范围 [最小值 最大值]

%% ========== 高斯滤波核设置 ==========
% 【优化说明】基于PS-OCT巩膜检测技术需求的参数调整：
% - 巩膜检测需要保持细节结构的同时抑制散斑噪声
% - 双折射测量要求在噪声抑制与空间分辨率间平衡
% - 考虑人眼巩膜组织的结构特征尺度
% - 基于DOPU图像质量反馈，增大滤波核以改善分层清晰度

% 小尺度高斯核 (用于细节保持和初步降噪) - 细节优先
params.filters.h1_size = [5 5];                % 高斯核1尺寸 (减小以增强局部纹理)
params.filters.h1_sigma = 1.5;                % 高斯核1标准差 (降低以避免过度平滑)
params.filters.h1 = fspecial('gaussian', params.filters.h1_size, params.filters.h1_sigma);

% 中尺度高斯核 (用于结构增强和背景平滑) - 背景优先
params.filters.h2_size = [19 19];               % 高斯核2尺寸 (显著增大以平滑大尺度背景)
params.filters.h2_sigma = 3;                  % 高斯核2标准差 (增大以抑制深层噪声)
params.filters.h2 = fspecial('gaussian', params.filters.h2_size, params.filters.h2_sigma);

% 【输出端自适应滤波参数】用于光轴和延迟结果的DOPU混合滤波策略
% 策略说明：
%   - DOPU >= 阈值（高质量组织）：使用固定h2滤波，保留细节
%   - DOPU < 阈值（低质量/深层）：使用自适应滤波，从h2核大小开始根据DOPU向上调整核大小
params.filters.enable_output_adaptive = 0;     % 1: 启用输出端混合滤波; 0: 使用传统固定h2滤波
params.filters.output_dopu_threshold = 0.35;    % DOPU阈值，区分高低质量区域（典型值0.3-0.5）
params.filters.kRL_output = 13;                % 自适应滤波核下限（通常设为h2核大小，作为低DOPU区域的起始核）
params.filters.kRU_output = 31;                % 自适应滤波核上限（DOPU=0时的最大核，用于低DOPU区域）
params.filters.adaptive_filter_bottom_depth = 80;  % DOPU自适应滤波仅对底部多少层生效（0表示全局生效，>0表示仅底部N层生效）

% 【底层相位延迟减小参数】用于深层区域的相位延迟降噪
% 在底层区域，对低DOPU像素降低相位延迟值以抑制噪声，提高深层成像质量
params.filters.enable_bottom_layer_phase_reduction = 0;  % 1: 启用底层相位延迟减小; 0: 禁用
params.filters.bottom_layer_depth = 120;        % 底层深度范围（从底部开始向上的层数）
params.filters.bottom_phase_reduction_ratio = 2;  % 相位延迟减小比例（0.5表示减小50%）
params.filters.bottom_dopu_threshold = 0.35;    % 底层DOPU阈值（小于此值的像素会被减小相位延迟）

% （已删除若干未使用的注释/备用参数，以保持配置简洁）

%% ========== 表面检测参数 ==========
% 【优化说明】针对眼部结构特点调整表面检测参数  
params.surface.surf_threshold = 0.20;          % 表面检测阈值(略降低以检测更多细微结构)
params.surface.median_window = 7;              % 中值滤波窗口大小(减小以保持边界清晰)
params.surface.smooth_window = 9;              % 平滑滤波窗口大小(减小以保持局部特征)

%% ========== 后巩膜边界处理参数 ==========
% 后巩膜边界文件路径（MAT文件）
% params.files.sclera_boundary_path = 'C:\yongxin.wang\Data\Process_Data\Disk\sclera_boundary.mat';  % 留空则不处理；设置完整路径则启用
% MAT文件中的变量名（可选，留空则自动探测）
params.files.sclera_boundary_var = '';   % 例如: 'boundary', 'layer', 'sclera_line' 等

%% ========== 日志与调试设置 ==========
% 全局日志等级: 0 = 静默 (默认), 1 = 简短, 2 = 详细
params.logging.verbose = 0;

%% ========== 验证参数合理性 ==========
% 【增强验证】检查所有关键参数的合理性
% 检查滤波核范围参数
if params.polarization.kRL_cfg1 >= params.polarization.kRU_cfg1
    warning('配置1滤波核范围设置有误: kRL_cfg1应小于kRU_cfg1');
    % 自动修正：交换值
    temp = params.polarization.kRL_cfg1;
    params.polarization.kRL_cfg1 = params.polarization.kRU_cfg1;
    params.polarization.kRU_cfg1 = temp;
    fprintf('已自动修正: kRL_cfg1=%d, kRU_cfg1=%d\n', params.polarization.kRL_cfg1, params.polarization.kRU_cfg1);
end

% 检查输出滤波核范围参数
if params.filters.kRL_output >= params.filters.kRU_output
    warning('输出滤波核范围设置有误: kRL_output应小于kRU_output');
    % 自动修正：交换值
    temp = params.filters.kRL_output;
    params.filters.kRL_output = params.filters.kRU_output;
    params.filters.kRU_output = temp;
    fprintf('已自动修正: kRL_output=%d, kRU_output=%d\n', params.filters.kRL_output, params.filters.kRU_output);
end

% （配置2 已弃用，相关检查移除）

% 调试输出：打印关键滤波参数
if params.logging.verbose >= 1
    fprintf('\n=== 关键滤波参数 ===\n');
    fprintf('kRL_cfg1 = %d, kRU_cfg1 = %d\n', params.polarization.kRL_cfg1, params.polarization.kRU_cfg1);
    fprintf('kRL_output = %d, kRU_output = %d\n', params.filters.kRL_output, params.filters.kRU_output);
    fprintf('=====================\n\n');
end

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
    fprintf('相位延迟范围 (PRRrg): [%.2f %.2f]\n', params.polarization.PRRrg(1), params.polarization.PRRrg(2));
    fprintf('展平并生成 En-face: %s\n', iif(params.processing.enable_flatten_enface, '启用', '禁用'));
    fprintf('DOPU相位抑制: %s\n', iif(params.polarization.enableDopuPhaseSupp, '启用', '禁用'));
    fprintf('输出端自适应滤波: %s', iif(params.filters.enable_output_adaptive, '启用', '禁用'));
    if params.filters.enable_output_adaptive
        fprintf(' (kR范围: %d~%d)\n', params.filters.kRL_output, params.filters.kRU_output);
    else
        fprintf('\n');
    end
    fprintf('底层相位延迟减小: %s', iif(params.filters.enable_bottom_layer_phase_reduction, '启用', '禁用'));
    if params.filters.enable_bottom_layer_phase_reduction
        fprintf(' (深度: %d层, 减小: %d%%, DOPU阈值: %.2f)\n', ...
            params.filters.bottom_layer_depth, ...
            round(params.filters.bottom_phase_reduction_ratio * 100), ...
            params.filters.bottom_dopu_threshold);
    else
        fprintf('\n');
    end
    fprintf('==========================\n\n');
end