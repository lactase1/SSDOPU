%% ====================================================================================
% 函数名: calculateSplitSpectrumDOPU_nonoverlap5
% 功能: 使用固定 5 段、无重叠的子频带计算分裂谱 DOPU
% 说明: 接口与 calculateSplitSpectrumDOPU 保持一致，便于在主脚本中互换使用
% 输入参数:
%   Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg
%   其中 nWin/windex/winL/winG 可由调用方传入，但本函数会忽略传入的 nWin/windex
%   并以 5 段无重叠策略重建自己的窗口参数
% 输出参数:
%   dopu_splitSpectrum - [nZcrop × nX]
%   dopu_ss - [nZcrop × nX × nr]
% 作者: 自动生成
% 日期: 2025
% ====================================================================================

function [dopu_splitSpectrum, dopu_ss] = calculateSplitSpectrumDOPU_nonoverlap5(Bd1, Bd2, params, SPL, nX, nr, ~, ~, ~, ~, czrg)

% 参数检查（与原函数行为相近）
[signalLength, nX_check, nr_check] = size(Bd1);
if signalLength ~= SPL || nX_check ~= nX || nr_check ~= nr
    error('输入信号尺寸与参数不匹配');
end

nZcrop = numel(czrg);

% 如果不开启分裂谱，返回全1（与原函数一致）
if ~isfield(params,'dopu') || ~isfield(params.dopu,'do_ssdopu') || ~params.dopu.do_ssdopu
    dopu_splitSpectrum = ones(nZcrop, nX);
    dopu_ss = ones(nZcrop, nX, nr);
    return;
end

% 固定 5 段，无重叠策略
fixed_nWin = 5;
% 为保证窗口对齐，我们按等分方式划分原始频谱长度 Blength=SPL
Blength = SPL;
winL = floor(Blength / fixed_nWin); % 每段长度，最后一段可能包含剩余
step = winL; % 无重叠
windex = 1:step:Blength;
% 若最后一个起点过多，裁剪到 fixed_nWin 个
if numel(windex) > fixed_nWin
    windex = windex(1:fixed_nWin);
end
% 构造窗函数（若 params 指定则使用，否则默认 tukey）
if isfield(params,'dopu') && isfield(params.dopu,'win_alpha')
    alpha = params.dopu.win_alpha;
else
    alpha = 0.25;
end

% 确保 winL 最少为 1
if winL < 1
    winL = 1;
end
winG = tukeywin(winL, alpha);

% 准备存储 FFT 结果
Bimg1 = zeros(SPL, nX, nr, numel(windex));
Bimg2 = zeros(SPL, nX, nr, numel(windex));

% 对每个重复和每个窗口进行 FFT（与原函数类似）
for iR = 1:nr
    for iL = 1:numel(windex)
        start_idx = windex(iL);
        end_idx = min(start_idx + winL - 1, signalLength);
        current_winL = end_idx - start_idx + 1;
        if current_winL ~= winL && length(winG) == winL
            current_winG = winG(1:current_winL);
        else
            current_winG = winG;
        end
        iBd1 = Bd1(start_idx:end_idx, :, iR) .* current_winG;
        iBd2 = Bd2(start_idx:end_idx, :, iR) .* current_winG;
        fft_result1 = fft(iBd1, SPL, 1);
        fft_result2 = fft(iBd2, SPL, 1);
        Bimg1(:, :, iR, iL) = fft_result1;
        Bimg2(:, :, iR, iL) = fft_result2;
    end
end

% 裁剪 Z
IMGs_ch1 = Bimg1(czrg, :, :, :);
IMGs_ch2 = Bimg2(czrg, :, :, :);

% 计算 Stokes
S0 = zeros(nZcrop, nX, nr, numel(windex));
S1 = zeros(nZcrop, nX, nr, numel(windex));
S2 = zeros(nZcrop, nX, nr, numel(windex));
S3 = zeros(nZcrop, nX, nr, numel(windex));
for iR = 1:nr
    for iL = 1:numel(windex)
        [S0(:,:,iR,iL), S1(:,:,iR,iL), S2(:,:,iR,iL), S3(:,:,iR,iL)] = ...
            cumulativeQUV(IMGs_ch1(:,:,iR,iL), IMGs_ch2(:,:,iR,iL));
    end
end

S1_avg = mean(S1./S0, 4);
S2_avg = mean(S2./S0, 4);
S3_avg = mean(S3./S0, 4);

% DOPU per repeat
dopu_ss = sqrt(S1_avg.^2 + S2_avg.^2 + S3_avg.^2);
% 平均重复
dopu_splitSpectrum = mean(dopu_ss, 3);

end
