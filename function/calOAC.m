%% ====================================================================================
% 函数名: calOAC
% 功能: 计算光学衰减系数(Optical Attenuation Coefficient)
% 输入参数:
%   linFrame - 线性OCT信号强度矩阵 [Z×X]，通常为|FFT(OCT信号)|²
% 输出参数:
%   OAC - 光学衰减系数矩阵 [Z×X]，用于组织边界检测和结构分析
% 算法原理:
%   - 基于Beer-Lambert定律: I(z) = I₀ exp(-2μz)
%   - μ为衰减系数，反映组织的散射和吸收特性
%   - 通过当前信号强度与剩余信号总和的比值估算局部衰减
% 应用: 主要用于组织表面检测、分层分析和结构图像增强
% ====================================================================================
function OAC = calOAC(linFrame)
    % 初始化输出
    OAC = zeros(size(linFrame));

    % 使用最后3层估算噪声底限
    last3 = linFrame(end-2:end, :);
    Nois_M = mean(last3(:));
    if Nois_M < 1
        return; % 信号过弱则直接返回
    end

    % 噪声标准差
    Nois_D = std(last3(:));

    % 去除噪声底限并添加最小阈值
    linFrame = linFrame - Nois_M + Nois_D;
    linFrame(linFrame < 0) = 1; % 避免对数或分母为零

    % 估算尾部残余光强 (最后6层)
    tail_signal = mean(linFrame(end-5:end, :), 1);
    tail_signal = medfilt1(tail_signal, 15);
    tail_signal = smooth(tail_signal, 15)';

    % 经验衰减参数 alpha
    alpha = 0.0086;

    % 计算每一深度层的 OAC
    Z = size(linFrame, 1);
    for z = 1:Z
        denom = 2 * alpha * sum(linFrame(z+1:Z, :), 1) + tail_signal;
        OAC(z, :) = linFrame(z, :) ./ denom;
    end
end
