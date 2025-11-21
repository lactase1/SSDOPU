%% vWinAvgFiltOpt: 动态加权高斯滤波函数
% 该函数根据输入权重动态选择高斯核大小，对输入数据进行局部加权高斯滤波。
% 适用于需要根据置信度调整滤波强度的场景，如偏振光学数据处理。
% 核心思想: 权重高时用小核 (精细滤波), 权重低时用大核 (强噪声抑制)。
%
% 输入参数:
%   inFrame - 输入数据矩阵 [nZ x nX]，待滤波的数据 (如 Stokes 参数)
%   inWeight - 权重矩阵 [nZ x nX]，每个像素的置信度 (0-1), 通常为 DOPU 值
%   kRL - 高斯核大小的下限 (整数), 最小核尺寸
%   kRU - 高斯核大小的上限 (整数), 最大核尺寸
%   Nsec - 分段数 (可选，默认 kRU - kRL + 1), 控制权重到核大小的映射精细度
%
% 输出参数:
%   outFrameWfilt - 滤波后的输出数据矩阵 [nZ x nX]
%
% 算法流程:
% 1. 根据权重范围生成一系列高斯核
% 2. 对每个像素，根据其权重选择合适的核
% 3. 使用选定核进行局部加权滤波
% 4. 处理边界通过数据扩展 (padding)
%
% 注意: 该版本的高斯核为矩形 [kRg(i), round(kRg(i)/2)], 且滤波时只用部分列,
%       这可能导致不对称滤波。建议使用 vWinAvgFiltOpt_2.m (方形核, 完整区域)。

function [outFrameWfilt] = vWinAvgFiltOpt(inFrame,inWeight,kRL,kRU,Nsec)
    % 处理可选参数 Nsec, 如果未提供则默认计算
    if nargin == 4, Nsec = kRU - kRL +1; end

    % 获取输入矩阵的尺寸: nZ (深度), nX (宽度)
    [nZ,nX] = size(inFrame);

    % 计算扩展边界的大小, 基于最大核尺寸 kRU
    % krs/kcs 为半径, 用于 padding
    krs = (kRU-1)/2; kcs = (kRU-1)/2;

    % 特殊情况: 如果核大小为0 (kRU=1), 直接返回简单加权结果, 无需滤波
    if krs == 0 && kcs ==0, outFrameWfilt = inFrame.*inWeight; return; end

    % 步骤2: 扩展输入数据和权重矩阵 (padding zeros)
    % 扩展大小为 krs/kcs, 以便边界像素也能正确滤波
    inframe = zeros(nZ+krs*2, nX+kcs*2);  % 扩展后的数据矩阵
    inweight = inframe;  % 扩展后的权重矩阵 (初始化为0)
    % 将原始数据复制到扩展矩阵的中心
    inframe(krs+1:end-krs, kcs+1:end-kcs) = inFrame;
    inweight(krs+1:end-krs, kcs+1:end-kcs) = inWeight;

    % 步骤3: 生成权重范围和对应的核大小
    % wRg: 权重范围 [0,1], 分成 Nsec 段
    wRg = linspace(0,1,Nsec);
    % kRg: 核大小范围 [kRU, kRL], 从大到小变化
    kRg = round(linspace(kRU,kRL,Nsec));

    % 生成一系列高斯核, 存储在 cell 数组 gaussHs 中
    % 核形状: 矩形 [kRg(i), round(kRg(i)/2)], 标准差 round(kRg(i)/2)
    for i = 1:1:Nsec
        gaussHs{i} = fspecial('gaussian', [kRg(i) round(kRg(i)/2)], round(kRg(i)/2));
    end

    % 步骤4: 遍历每个像素，进行动态加权滤波
    for iz = 1:nZ  % 深度循环
        for ix = 1:nX  % 宽度循环
            % 根据当前像素的权重 inWeight(iz,ix) 选择合适的核索引
            % wInds: 布尔数组, 标记权重 >= 当前值的索引
            wInds = wRg >= inWeight(iz,ix);
            pInds = find(wInds);  % 找到满足条件的索引
            fInd = pInds(1);  % 选择第一个 (权重高时选小索引, 小核)
            fkrg = kRg(fInd);  % 对应的核大小

            % 计算核的半径范围, 处理奇偶性
            if mod(fkrg,2)  % 奇数核
                ks1 = (fkrg-1)/2; ks2 = ks1;
            else  % 偶数核
                ks1 = fkrg/2-1;
                ks2 = ks1+1;
            end

            % 定义局部区域的行列索引 (相对于扩展矩阵)
            cols = ix + kcs - ks1 : ix + kcs + ks2;  % 列索引
            rows = iz + krs - ks1 : iz + krs + ks2;  % 行索引

            % 提取局部数据 (注意: 只取部分列 cols(1:round(fkrg/2)), 可能不对称!)
            inA = inframe(rows, cols(1:round(fkrg/2)));

            % 使用选定的高斯核进行加权滤波: 数据 * 核, 求和
            outFrameWfilt(iz,ix) = sum(inA .* gaussHs{fInd}, 'all');
        end
    end
end