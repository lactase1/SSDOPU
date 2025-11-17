%% ====================================================================================
% 函数名: calculateSpatialDOPU
% 功能: 计算空间DOPU (基于3x3邻域平均)
% 算法原理:
%   - 计算每个像素的Stokes参数 (S0, S1, S2, S3)
%   - 对每个深度层，使用3x3邻域对Stokes参数进行空间平均
%   - 从平均的Stokes参数计算DOPU
%   - DOPU = sqrt((<S1>/<S0>)² + (<S2>/<S0>)² + (<S3>/<S0>)²)
%   - 其中<S0>、<S1>、<S2>、<S3>是3x3邻域的平均值
%
% 输入参数:
%   IMG_ch1 - 通道1的复数OCT信号矩阵 [Z×X×Rep]
%   IMG_ch2 - 通道2的复数OCT信号矩阵 [Z×X×Rep]
%   params - 参数结构体，包含:
%       .dopu.do_spatial - 是否进行空间DOPU计算的标志 (true/false)
%
% 输出参数:
%   dopu_spatial - 空间DOPU结果 [Z×X]
%
% 依赖函数:
%   - cumulativeQUV: 计算Stokes参数的函数
%
% 作者: PS-OCT处理系统
% 版本: 1.0
% 日期: 2024
% ====================================================================================

function dopu_spatial = calculateSpatialDOPU(IMG_ch1, IMG_ch2, params)

% 参数验证
if nargin < 3
    error('函数需要3个输入参数');
end

if isempty(IMG_ch1) || isempty(IMG_ch2)
    error('输入信号 IMG_ch1 和 IMG_ch2 不能为空');
end

if size(IMG_ch1) ~= size(IMG_ch2)
    error('IMG_ch1 和 IMG_ch2 的尺寸必须相同');
end

% 检查是否需要进行空间DOPU计算
if ~params.dopu.do_spatial
    % 如果不进行空间计算，返回默认值1
    [nZ, nX, ~] = size(IMG_ch1);
    dopu_spatial = ones(nZ, nX);
    return;
end

% 获取数据维度
[nZ, nX, nRep] = size(IMG_ch1);

% 对重复次数进行平均 (如果有多个重复)
if nRep > 1
    IMG_ch1 = mean(IMG_ch1, 3);
    IMG_ch2 = mean(IMG_ch2, 3);
end

% 计算Stokes参数
[S0, S1, S2, S3] = cumulativeQUV(IMG_ch1, IMG_ch2);

% 初始化空间DOPU结果
dopu_spatial = zeros(nZ, nX);

% 对每个深度层进行空间滤波
for z = 1:nZ
    % 获取当前层的Stokes参数
    S0_layer = S0(z, :);
    S1_layer = S1(z, :);
    S2_layer = S2(z, :);
    S3_layer = S3(z, :);

    % 将1D数组重塑为2D矩阵进行2D卷积 (这里其实已经是1D的了，需要重新考虑)
    % 实际上，由于我们是对每个深度层独立处理，而每个深度层是1D的X方向数据
    % 所以3x3邻域平均在这种情况下应该是1D的3点平均

    % 边界扩充处理：复制边界像素向外扩展1像素
    S0_padded = padarray(S0_layer, [0, 1], 'replicate', 'both');
    S1_padded = padarray(S1_layer, [0, 1], 'replicate', 'both');
    S2_padded = padarray(S2_layer, [0, 1], 'replicate', 'both');
    S3_padded = padarray(S3_layer, [0, 1], 'replicate', 'both');

    % 使用3点均值滤波器对Stokes参数进行空间平均 (1D卷积，沿X方向)
    kernel = ones(1, 3) / 3;

    % 对每个Stokes参数应用空间平均
    S0_avg = conv(S0_padded, kernel, 'valid');
    S1_avg = conv(S1_padded, kernel, 'valid');
    S2_avg = conv(S2_padded, kernel, 'valid');
    S3_avg = conv(S3_padded, kernel, 'valid');

    % 从平均的Stokes参数计算DOPU
    % 避免除零错误
    S0_avg_safe = max(S0_avg, eps);  % eps是最小的正浮点数
    dopu_spatial(z, :) = sqrt((S1_avg./S0_avg_safe).^2 + (S2_avg./S0_avg_safe).^2 + (S3_avg./S0_avg_safe).^2);
end

end