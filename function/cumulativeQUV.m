%% ====================================================================================
% 函数名: cumulativeQUV
% 功能: 从两通道复数OCT图像计算Stokes参数(S0, S1, S2, S3)
% 算法原理:
%   - Stokes参数是描述偏振光的完整参数集
%   - S0: 总光强度 |E1|² + |E2|²
%   - S1: 水平/垂直偏振差 |E1|² - |E2|²
%   - S2: ±45°偏振差 2|E1||E2|cos(θ)
%   - S3: 左/右圆偏振差 2|E1||E2|sin(-θ)
%   - θ为两通道信号的相位差
%
% 输入参数:
%   IMG_ch1 - 通道1的复数OCT信号矩阵 [Z×X或Z×X×Rep]
%   IMG_ch2 - 通道2的复数OCT信号矩阵 [Z×X或Z×X×Rep]
%
% 输出参数:
%   S0 - Stokes参数S0 (总光强度)
%   S1 - Stokes参数S1 (水平/垂直偏振差)
%   S2 - Stokes参数S2 (±45°偏振差)
%   S3 - Stokes参数S3 (左/右圆偏振差)
%
% 调用示例:
%   [S0, S1, S2, S3] = cumulativeQUV(IMG_ch1, IMG_ch2);
%
% 应用: 用于DOPU(偏振均匀度)和各种偏振参数计算
% ====================================================================================

function [S0, S1, S2, S3] = cumulativeQUV(IMG_ch1, IMG_ch2)
% 计算两通道信号的相位差
axis = angle(IMG_ch2 .* conj(IMG_ch1));

% 计算四个Stokes参数
S0 = abs(IMG_ch1).^2 + abs(IMG_ch2).^2;                    % 总光强度
S1 = abs(IMG_ch1).^2 - abs(IMG_ch2).^2;                    % 水平-垂直偏振分量差
S2 = 2 .* abs(IMG_ch1) .* abs(IMG_ch2) .* cos(axis);       % +45°与-45°偏振分量差
S3 = 2 .* abs(IMG_ch1) .* abs(IMG_ch2) .* sin(-axis);      % 右旋与左旋圆偏振分量差
end