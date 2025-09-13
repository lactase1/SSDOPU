function [outFrameWfilt] = vWinAvgFiltOpt_2_1(inFrame, inWeight, kRL, kRU, Nsec)
%% VWINAVGFILTOPT_2_1 自适应可变窗口高斯滤波nction [outFrameWfilt] = vWinAvgFiltOpt_2_1(inFrame, inWeight, kRL, kRU, Nsec)
%   [outFrameWfilt] = vWinAvgFiltOpt_2_1(inFrame, inWeight, kRL, kRU, Nsec)
%
%   本函数实现基于置信度的自适应可变窗口高斯滤波，可用于PSOCT图像处理
%   核心思想：高置信度区域使用小窗口保留细节，低置信度区域使用大窗口增强平滑
%
% 输入参数:
%   inFrame  - 需要滤波的二维数据矩阵 [nZ×nX]
%   inWeight - 置信度矩阵（如DOPU值）[nZ×nX]，范围0-1
%              值越大表示置信度越高，会使用越小的滤波窗口
%   kRL      - 高斯核最小尺寸(下限)，用于高置信度区域，推荐值2-5
%   kRU      - 高斯核最大尺寸(上限)，用于低置信度区域，推荐值7-21
%   Nsec     - 窗口大小划分的区间数，默认为kRU-kRL+1
%
% 输出参数:
%   outFrameWfilt - 滤波后的二维数据矩阵 [nZ×nX]
%
% 参数调整建议:
%   - 对于需要保留精细结构的区域，使用较小的kRL (2-3)
%   - 对于需要强平滑的区域，使用较大的kRU (15-21)
%   - 对于一般平衡需求，推荐kRL=3, kRU=9
%
% 示例:
%   outFrame = vWinAvgFiltOpt_2_1(Q, DOPU, 3, 9);  % 使用DOPU作为置信度滤波Q矩阵
%
% 相关函数:
%   fspecial, imfilter

    % 设置默认Nsec值（如果未提供）
    if nargin == 4
        Nsec = kRU - kRL + 1;
    end
    
    %% 1. 初始化和边界处理
    % 获取输入数据的尺寸
    [nZ, nX] = size(inFrame);
    
    % 计算最大窗口半径（用于边界填充）
    krs = (kRU-1)/2; 
    kcs = (kRU-1)/2;
    
    % 特殊情况处理：如果窗口大小为1x1，则直接返回结果
    if krs == 0 && kcs == 0
        outFrameWfilt = inFrame.*inWeight; 
        return;
    end
    
    % 数据扩展：创建填充边界的扩展矩阵
    inframe = zeros(nZ+krs*2, nX+kcs*2); 
    
    % 将原始数据放入扩展矩阵的中心区域
    inframe(krs+1:end-krs, kcs+1:end-kcs) = inFrame;
    
    % 同样扩展权重矩阵（虽然代码中没有使用）
    inweight = zeros(nZ+krs*2, nX+kcs*2);
    inweight(krs+1:end-krs, kcs+1:end-kcs) = inWeight;
    
    %% 2. 根据权重阈值生成高斯核集合
    % 生成权重阈值序列，均匀分布在0到1之间
    wRg = linspace(0, 1, Nsec);
    
    % 生成对应的窗口大小序列，从大到小（kRU到kRL）
    kRg = round(linspace(kRU, kRL, Nsec));
    
    % 预生成所有可能使用的高斯核
    gaussHs = cell(1, Nsec);
    for i = 1:Nsec
        % 创建椭圆形高斯核: [高度 宽度为高度的一半]，标准差为尺寸的一半
        gaussHs{i} = fspecial('gaussian', [kRg(i) round(kRg(i)/2)], round(kRg(i)/2));
    end
    
    %% 3. 对每个像素进行自适应滤波
    % 预分配输出矩阵
    outFrameWfilt = zeros(nZ, nX);
    
    % 对每个像素循环处理
    for iz = 1:nZ
        for ix = 1:nX
            % 根据当前像素的置信度值选择合适的窗口大小
            % 找出所有大于等于当前置信度的阈值索引
            wInds = wRg >= inWeight(iz, ix);
            pInds = find(wInds);
            
            % 选择第一个符合条件的窗口大小（注：越小的索引对应越大的窗口）
            fInd = pInds(1); 
            fkrg = kRg(fInd);
            
            % 计算窗口半径（处理奇偶不同情况）
            if mod(fkrg, 2)  % 奇数大小窗口
                ks1 = (fkrg-1)/2; 
                ks2 = ks1;
            else  % 偶数大小窗口
                ks1 = fkrg/2-1; 
                ks2 = ks1+1;
            end
            
            % 计算当前像素邻域的索引范围
            cols = ix+kcs-ks1:ix+kcs+ks2;
            rows = iz+krs-ks1:iz+krs+ks2;
            
            % 提取当前窗口区域的数据（调整为高斯核的大小）
            inA = inframe(rows, cols(1:round(fkrg/2)));
            
            % 应用高斯核进行加权平均，得到滤波后的值
            outFrameWfilt(iz, ix) = sum(inA.*gaussHs{fInd}, 'all');
        end
    end
end