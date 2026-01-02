function V_out = filterVectorFieldAdaptive(V_in, weightMap, kRL, kRU, h2, dopuThreshold)
% filterVectorFieldAdaptive - 基于DOPU阈值的混合向量场滤波
%
% 功能说明:
%   对三维向量场（如光轴LA）进行基于DOPU阈值的分区域滤波策略：
%   - DOPU >= 阈值（高质量区域）：使用固定h2高斯滤波，保留组织细节
%   - DOPU < 阈值（低质量区域）：使用DOPU自适应滤波，从h2核大小开始根据DOPU值向上调整核大小以抑制噪声
%   - 滤波后重新归一化，确保向量单位长度
%
% 输入参数:
%   V_in         - 输入向量场 [nZ, nX, 3]，第三维为向量的三个分量
%   weightMap    - DOPU权重图 [nZ, nX]，范围[0,1]，值越大表示信号质量越好
%   kRL          - 自适应滤波核半径下限（对应DOPU接近阈值时的核大小，通常设为h2核大小）
%   kRU          - 自适应滤波核半径上限（对应DOPU=0时的最大核）
%   h2           - 固定高斯滤波核，用于高DOPU区域 [可选，如不提供则kRL将作为固定核大小]
%   dopuThreshold - DOPU阈值，用于区分高低质量区域 [可选，默认0.4]
%
% 输出参数:
%   V_out        - 输出向量场 [nZ, nX, 3]，归一化后的混合滤波结果
%
% 算法逻辑:
%   1. 创建DOPU阈值掩码
%   2. 高DOPU区域：对三个分量分别使用h2固定高斯滤波
%   3. 低DOPU区域：对三个分量使用DOPU自适应滤波（vWinAvgFiltOpt）
%   4. 合并两部分结果并重新归一化
%
% 使用示例:
%   cumLA_filtered = filterVectorFieldAdaptive(cumLA_raw, DopuF, 13, 21, h2, 0.4);
%
% 创建时间: 2025年12月31日
% 修改时间: 2026年1月1日
% 作者: Yongxin Wang
% 关联功能: DOPU混合滤波策略（固定+自适应）

%% 参数验证
if nargin < 4
    error('filterVectorFieldAdaptive: 至少需要4个输入参数 (V_in, weightMap, kRL, kRU)');
end

[nZ, nX, nDim] = size(V_in);
if nDim ~= 3
    error('filterVectorFieldAdaptive: V_in的第三维必须为3（向量场）');
end

if any(size(weightMap) ~= [nZ, nX])
    error('filterVectorFieldAdaptive: weightMap尺寸必须与V_in的前两维匹配');
end

if kRL >= kRU
    warning('filterVectorFieldAdaptive: kRL应小于kRU，自动修正为 kRU=kRL+1');
    kRU = kRL + 1;
end

% 如果没有提供h2，则使用kRL大小创建高斯核
if nargin < 5 || isempty(h2)
    h2 = fspecial('gaussian', [kRL kRL], max(1, round(kRL/2)));
end

% 如果没有提供DOPU阈值，则使用默认值0.4
if nargin < 6 || isempty(dopuThreshold)
    dopuThreshold = 0.4;
end

%% DOPU阈值掩码
% 根据配置的阈值区分高低质量区域
maskHighQuality = (weightMap >= dopuThreshold);  % 高DOPU区域（使用固定h2滤波）
maskLowQuality = (weightMap > 0) & (weightMap < dopuThreshold);  % 低DOPU区域（使用自适应滤波）

%% 分区域滤波
% 初始化输出
V1 = zeros(nZ, nX);
V2 = zeros(nZ, nX);
V3 = zeros(nZ, nX);

% === 区域1: 高DOPU区域（>= 0.4）使用固定h2滤波 ===
if any(maskHighQuality(:))
    V1_high = imfilter(V_in(:,:,1), h2, 'replicate');
    V2_high = imfilter(V_in(:,:,2), h2, 'replicate');
    V3_high = imfilter(V_in(:,:,3), h2, 'replicate');
    
    V1(maskHighQuality) = V1_high(maskHighQuality);
    V2(maskHighQuality) = V2_high(maskHighQuality);
    V3(maskHighQuality) = V3_high(maskHighQuality);
end

% === 区域2: 低DOPU区域（< 0.4）使用自适应DOPU滤波 ===
if any(maskLowQuality(:))
    % 从kRL（h2核大小）开始，根据DOPU值向上调整到kRU
    % DOPU越低，核越大，滤波越强
    V1_low = vWinAvgFiltOpt(V_in(:,:,1), weightMap, kRL, kRU);
    V2_low = vWinAvgFiltOpt(V_in(:,:,2), weightMap, kRL, kRU);
    V3_low = vWinAvgFiltOpt(V_in(:,:,3), weightMap, kRL, kRU);
    
    V1(maskLowQuality) = V1_low(maskLowQuality);
    V2(maskLowQuality) = V2_low(maskLowQuality);
    V3(maskLowQuality) = V3_low(maskLowQuality);
end

% === 区域3: DOPU=0的无效区域保持原值 ===
maskInvalid = (weightMap == 0);
if any(maskInvalid(:))
    V1(maskInvalid) = V_in(maskInvalid, 1);
    V2(maskInvalid) = V_in(maskInvalid, 2);
    V3(maskInvalid) = V_in(maskInvalid, 3);
end

%% 向量归一化
% 【关键步骤】滤波后向量模长改变，必须重新归一化
% 计算向量的模长（欧几里得范数）
normV = sqrt(V1.^2 + V2.^2 + V3.^2);

% 防止除以零：将模长为0的位置设为1（这些位置会在后续被NaN处理掉）
normV(normV < 1e-9) = 1;

% 归一化：确保输出向量的模长为1
V_out = cat(3, V1./normV, V2./normV, V3./normV);

%% 数值稳定性处理
% 将NaN值替换为0，避免后续计算出现问题
V_out(isnan(V_out)) = 0;
V_out(isinf(V_out)) = 0;

end
