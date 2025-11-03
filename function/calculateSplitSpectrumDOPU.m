%% ====================================================================================
% 函数名: calculateSplitSpectrumDOPU
% 功能: 计算分裂谱偏振均匀度(Degree of Polarization Uniformity)
% 算法原理:
%   - 将OCT信号频谱分割成多个子频带
%   - 对每个子频带分别计算偏振参数(Stokes参数)
%   - 通过频谱平均提高信噪比和空间分辨率
%   - 计算DOPU = sqrt(<S1/S0>² + <S2/S0>² + <S3/S0>²)
%
% 输入参数:
%   Bd1 - 通道1的OCT信号矩阵 [信号长度×X像素×重复次数]
%   Bd2 - 通道2的OCT信号矩阵 [信号长度×X像素×重复次数]
%   params - 参数结构体，包含:
%       .dopu.do_ssdopu - 是否进行分裂谱计算的标志 (true/false)
%   SPL - 信号长度 (FFT点数)
%   nX - X方向像素数
%   nr - 重复次数 (repetitions)
%   nWin - 分裂谱窗口数 (默认9)
%   windex - 窗口起始索引数组 [1×nWin]
%   winL - 单个窗口长度
%   winG - 高斯窗口函数 [winL×1]
%   czrg - Z方向裁剪范围 (如: 1:320)
%
% 输出参数:
%   dopu_splitSpectrum - 分裂谱DOPU结果 [Z像素×X像素]
%   dopu_ss - 各重复次数的分裂谱DOPU [Z像素×X像素×重复次数]
%
% 调用示例:
%   % 基本调用
%   [dopu_result, dopu_all] = calculateSplitSpectrumDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);
%
%   % 如果不进行分裂谱计算
%   params.dopu.do_ssdopu = false;
%   [dopu_result, ~] = calculateSplitSpectrumDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);
%   % 返回 dopu_result = 1 (默认值)
%
% 依赖函数:
%   - cumulativeQUV: 计算Stokes参数的函数
%
% 作者: PS-OCT处理系统
% 版本: 1.0
% 日期: 2024
% ====================================================================================

function [dopu_splitSpectrum, dopu_ss] = calculateSplitSpectrumDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg)

% 参数验证
if nargin < 11
    error('函数需要11个输入参数');
end

if isempty(Bd1) || isempty(Bd2)
    error('输入信号 Bd1 和 Bd2 不能为空');
end

if size(Bd1) ~= size(Bd2)
    error('Bd1 和 Bd2 的尺寸必须相同');
end

% 获取矩阵尺寸
[signalLength, nX_check, nr_check] = size(Bd1);

if signalLength ~= SPL || nX_check ~= nX || nr_check ~= nr
    error('输入信号尺寸与参数不匹配');
end

% 计算Z方向裁剪后的长度
nZcrop = numel(czrg);

% 初始化输出变量
dopu_splitSpectrum = [];
dopu_ss = [];

% 控制打印冗余信息：优先使用 params.dopu.verbose，其次回退到 params.processing.show_img，默认 false
if isfield(params,'dopu') && isfield(params.dopu,'verbose')
    verbose = logical(params.dopu.verbose);
elseif isfield(params,'processing') && isfield(params.processing,'show_img')
    verbose = logical(params.processing.show_img);
else
    verbose = false;
end

% 检查是否需要进行分裂谱计算
if ~params.dopu.do_ssdopu
    % 如果不进行分裂谱计算，返回默认值1
    dopu_splitSpectrum = ones(nZcrop, nX);
    dopu_ss = ones(nZcrop, nX, nr);
    return;
end

%% 核心分裂谱DOPU计算

% 步骤1: 创建数组存储分裂谱复数FFT结果 [Z×X×重复次数×窗口数]
if verbose
    fprintf('  分裂谱DOPU: 窗口数=%d, 重复=%d, X=%d, 裁剪Z=%d\n', nWin, nr, nX, nZcrop);
end

Bimg1 = zeros(SPL, nX, nr, nWin);  % 通道1的FFT结果
Bimg2 = zeros(SPL, nX, nr, nWin);  % 通道2的FFT结果

% 步骤2: 对每个重复和每个窗口进行FFT
for iR = 1:nr
    for iL = 1:nWin
        % 提取当前窗口的数据片段并应用高斯窗口
        start_idx = windex(iL);
        end_idx = min(start_idx + winL - 1, signalLength);  % 防止越界

        % 确保窗口长度正确
        current_winL = end_idx - start_idx + 1;
        if current_winL ~= winL && length(winG) == winL
            % 如果窗口长度不匹配，调整高斯窗口
            current_winG = winG(1:current_winL);
        else
            current_winG = winG;
        end

        % 应用窗口函数
        iBd1 = Bd1(start_idx:end_idx, :, iR) .* current_winG;
        iBd2 = Bd2(start_idx:end_idx, :, iR) .* current_winG;

        % 执行FFT并存储结果
        fft_result1 = fft(iBd1, SPL, 1);  % [SPL × nX]
        fft_result2 = fft(iBd2, SPL, 1);  % [SPL × nX]

        % 存储FFT结果
        Bimg1(:, :, iR, iL) = fft_result1;
        Bimg2(:, :, iR, iL) = fft_result2;
    end
end

% 步骤3: 裁剪Z方向范围
IMGs_ch1 = Bimg1(czrg, :, :, :);  % [nZcrop×nX×nr×nWin]
IMGs_ch2 = Bimg2(czrg, :, :, :);

% 步骤4: 计算每个窗口和重复的Stokes参数

S0 = zeros(nZcrop, nX, nr, nWin);
S1 = zeros(nZcrop, nX, nr, nWin);
S2 = zeros(nZcrop, nX, nr, nWin);
S3 = zeros(nZcrop, nX, nr, nWin);

for iR = 1:nr
    for iL = 1:nWin
        % 计算当前窗口和重复的Stokes参数
        [S0(:,:,iR,iL), S1(:,:,iR,iL), S2(:,:,iR,iL), S3(:,:,iR,iL)] = ...
            cumulativeQUV(IMGs_ch1(:,:,iR,iL), IMGs_ch2(:,:,iR,iL));
    end
end

% 步骤5: 计算分裂谱DOPU

% 对所有窗口进行平均，然后计算DOPU
% S1./S0, S2./S0, S3./S0 分别在第4维度(窗口)上平均
S1_avg = mean(S1./S0, 4);  % [nZcrop×nX×nr]
S2_avg = mean(S2./S0, 4);
S3_avg = mean(S3./S0, 4);

% 计算DOPU: sqrt(<S1/S0>² + <S2/S0>² + <S3/S0>²)
dopu_ss = sqrt(S1_avg.^2 + S2_avg.^2 + S3_avg.^2);  % [nZcrop×nX×nr]

% 对所有重复进行平均，得到最终结果
dopu_splitSpectrum = mean(dopu_ss, 3);  % [nZcrop×nX]

if verbose
    fprintf('  分裂谱DOPU计算完成\n');
end

end