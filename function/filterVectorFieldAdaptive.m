function V_out = filterVectorFieldAdaptive(V_in, weightMap, kRL, kRU)
% filterVectorFieldAdaptive - 基于DOPU权重的自适应向量场滤波
%
% 功能说明:
%   对三维向量场（如光轴LA）进行基于DOPU的自适应平滑滤波。
%   - 高DOPU区域（组织）：小核滤波，保留细节
%   - 低DOPU区域（噪声/背景）：大核滤波，强平滑抑制噪声
%   - 滤波后重新归一化，确保向量单位长度
%
% 输入参数:
%   V_in      - 输入向量场 [nZ, nX, 3]，第三维为向量的三个分量
%   weightMap - DOPU权重图 [nZ, nX]，范围[0,1]，值越大表示信号质量越好
%   kRL       - 滤波核半径下限（高DOPU区域使用小核）
%   kRU       - 滤波核半径上限（低DOPU区域使用大核）
%
% 输出参数:
%   V_out     - 输出向量场 [nZ, nX, 3]，归一化后的自适应滤波结果
%
% 算法逻辑:
%   1. 分别对向量的三个分量进行DOPU自适应滤波（调用vWinAvgFiltOpt）
%   2. 滤波后向量模长不再为1，需要重新归一化
%   3. 处理零向量和NaN值，确保数值稳定性
%
% 使用示例:
%   cumLA_filtered = filterVectorFieldAdaptive(cumLA_raw, DopuF, 3, 21);
%
% 创建时间: 2025年12月31日
% 作者: Yongxin Wang
% 关联功能: DOPU自适应输出滤波（光轴和延迟结果后处理）

%% 参数验证
if nargin < 4
    error('filterVectorFieldAdaptive: 需要4个输入参数 (V_in, weightMap, kRL, kRU)');
end

[nZ, nX, nDim] = size(V_in);
if nDim ~= 3
    error('filterVectorFieldAdaptive: V_in的第三维必须为3（向量场）');
end

if any(size(weightMap) ~= [nZ, nX])
    error('filterVectorFieldAdaptive: weightMap尺寸必须与V_in的前两维匹配');
end

if kRL >= kRU
    warning('filterVectorFieldAdaptive: kRL应小于kRU，使用固定核kRL进行滤波');
    kRU = kRL + 1;
end

%% 分量自适应滤波
% 对向量场的三个分量分别进行DOPU加权自适应滤波
% 【关键】这里复用已有的vWinAvgFiltOpt函数，实现高效的自适应滤波
V1 = vWinAvgFiltOpt(V_in(:,:,1), weightMap, kRL, kRU);
V2 = vWinAvgFiltOpt(V_in(:,:,2), weightMap, kRL, kRU);
V3 = vWinAvgFiltOpt(V_in(:,:,3), weightMap, kRL, kRU);

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
