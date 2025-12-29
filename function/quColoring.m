%% ====================================================================================
% 函数名: quColoring (QU彩色编码)
% 功能: 将双折射的Q、U分量转换为HSV彩色编码图像
% 输入参数:
%   fLAs - 双折射数据矩阵 [nZ×nX×3×nY]，第3维为[?, Q, U]分量
%   cmapshift - 颜色映射旋转角度(可选参数)，用于调整色彩显示
% 输出参数:
%   hsvLA - HSV彩色编码图像 [nZ×nX×3×nY]，第3维为[H,S,V]分量
% 算法原理:
%   - 使用Q、U分量计算角度: θ = atan2(U, Q)
%   - 角度映射到色调(Hue): θ ∈ [-π, π] → H ∈ [0, 512]
%   - 使用HSV颜色空间表示双折射方向和强度
%   - 支持颜色映射的循环移位以优化显示效果
% 应用场景:
%   - 双折射方向的可视化
%   - 胶原纤维取向分析
%   - 病理组织的偏振特征显示
% 颜色含义:
%   - 色调(H): 双折射轴的方向角
%   - 饱和度(S): 通常设为最大值
%   - 亮度(V): 可用于表示双折射强度
% ====================================================================================
function [hsvLA] = quColoring(fLAs, cmapshift)
    % 创建HSV颜色映射表(512色)
    cusmapRg = hsv(512);

    % 如果提供了颜色映射移位参数，则循环移位颜色表
    if nargin > 1
        cusmapRg = circshift(cusmapRg, cmapshift);
    end

    % 定义角度范围[-π, π]对应的索引范围[1, 512]
    thetaRg = linspace(-pi, pi, 256*2);

    % 计算双折射轴角度: θ = atan2(U, Q)
    % fLAs(:,:,2) = U分量, fLAs(:,:,1) = Q分量
    thetas = atan2(fLAs(:,:,2), fLAs(:,:,1));

    % 使用多项式拟合建立角度到索引的映射关系
    thpf = polyfit(thetaRg, 1:256*2, 1);
    thetasInds = round(polyval(thpf, thetas));
    
    % 根据索引从颜色映射表中提取颜色
    colorInds = cusmapRg(thetasInds, :);
    
    % 重塑为与输入相同的维度结构
    hsvLA = reshape(colorInds, [size(thetasInds), 3]);
end
