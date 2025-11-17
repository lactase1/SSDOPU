%% ====================================================================================
% 函数名: calculateCombinedDOPU
% 功能: 计算分裂谱+空间组合DOPU
% 算法原理:
%   1. 先进行分裂谱DOPU计算，得到每个像素的分裂谱DOPU
%   2. 然后对分裂谱DOPU结果进行空间滤波(3x3邻域平均)
%   3. 得到最终的组合DOPU结果
%
% 输入参数:
%   Bd1 - 通道1的OCT信号矩阵 [信号长度×X像素×重复次数]
%   Bd2 - 通道2的OCT信号矩阵 [信号长度×X像素×重复次数]
%   params - 参数结构体，包含:
%       .dopu.do_combined - 是否进行组合DOPU计算的标志 (true/false)
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
%   dopu_combined - 组合DOPU结果 [Z像素×X像素]
%
% 依赖函数:
%   - calculateSplitSpectrumDOPU: 分裂谱DOPU计算函数
%
% 作者: PS-OCT处理系统
% 版本: 1.0
% 日期: 2024
% ====================================================================================

function dopu_combined = calculateCombinedDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg)

% 参数验证
if nargin < 11
    error('函数需要11个输入参数');
end

% 检查是否需要进行组合DOPU计算
if ~params.dopu.do_combined
    % 如果不进行组合计算，返回默认值1
    nZcrop = numel(czrg);
    dopu_combined = ones(nZcrop, nX);
    return;
end

% 步骤1: 计算分裂谱DOPU
[dopu_split, ~] = calculateSplitSpectrumDOPU(Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);

% 步骤2: 对分裂谱DOPU进行空间滤波
% 初始化组合DOPU结果
[nZ, nX] = size(dopu_split);
dopu_combined = zeros(nZ, nX);

% 对每个深度层进行空间滤波
for z = 1:nZ
    % 获取当前层的分裂谱DOPU
    dopu_layer = dopu_split(z, :);

    % 边界扩充处理：复制边界像素向外扩展1像素
    dopu_padded = padarray(dopu_layer, [0, 1], 'replicate', 'both');

    % 使用3x3均值滤波器进行空间平均 (实际上是1D的3点平均)
    % 创建3点均值核
    kernel = ones(1, 3) / 3;

    % 应用卷积进行空间平均
    dopu_filtered = conv(dopu_padded, kernel, 'valid');

    % 存储结果
    dopu_combined(z, :) = dopu_filtered;
end

end