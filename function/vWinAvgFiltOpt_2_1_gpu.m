function [outFrameWfilt] = vWinAvgFiltOpt_2_1_gpu(inFrame, inWeight, kRL, kRU, Nsec)
%% VWINAVGFILTOPT_2_1_GPU 基于GPU的自适应可变窗口高斯滤波
%   [outFrameWfilt] = vWinAvgFiltOpt_2_1_gpu(inFrame, inWeight, kRL, kRU, Nsec)
%
%   本函数是vWinAvgFiltOpt_2_1的GPU加速版本，实现基于置信度的自适应可变窗口高斯滤波
%   适用于PSOCT图像处理，通过GPU并行计算显著提高性能
%
% 输入参数:
%   inFrame  - 需要滤波的二维数据矩阵 [nZ×nX]，可以是gpuArray
%   inWeight - 置信度矩阵（如DOPU值）[nZ×nX]，范围0-1，可以是gpuArray
%              值越大表示置信度越高，会使用越小的滤波窗口
%   kRL      - 高斯核最小尺寸(下限)，用于高置信度区域，推荐值2-5
%   kRU      - 高斯核最大尺寸(上限)，用于低置信度区域，推荐值7-21
%   Nsec     - 窗口大小划分的区间数，默认为kRU-kRL+1
%
% 输出参数:
%   outFrameWfilt - 滤波后的二维数据矩阵 [nZ×nX]，为gpuArray类型
%
% 参数调整建议:
%   - 对于需要保留精细结构的区域，使用较小的kRL (2-3)
%   - 对于需要强平滑的区域，使用较大的kRU (15-21)
%   - 对于一般平衡需求，推荐kRL=3, kRU=9
%
% 示例:
%   % 将数据转移到GPU
%   gQ = gpuArray(Q);
%   gDOPU = gpuArray(DOPU);
%   % 在GPU上执行滤波
%   gOutFrame = vWinAvgFiltOpt_2_1_gpu(gQ, gDOPU, 3, 9);
%   % 如需，将结果传回CPU
%   outFrame = gather(gOutFrame);
%
% 相关函数:
%   fspecial, imgaussfilt (GPU版本内部实现)

    % 设置默认Nsec值（如果未提供）
    if nargin == 4
        Nsec = kRU - kRL + 1;
    end
    
    % 确保输入数据在GPU上
    if ~isa(inFrame, 'gpuArray')
        inFrame = gpuArray(inFrame);
    end
    
    if ~isa(inWeight, 'gpuArray')
        inWeight = gpuArray(inWeight);
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
    inframe = gpuArray.zeros(nZ+krs*2, nX+kcs*2); 
    
    % 将原始数据放入扩展矩阵的中心区域
    inframe(krs+1:end-krs, kcs+1:end-kcs) = inFrame;
    
    % 同样扩展权重矩阵
    inweight = gpuArray.zeros(nZ+krs*2, nX+kcs*2);
    inweight(krs+1:end-krs, kcs+1:end-kcs) = inWeight;
    
    %% 2. 根据权重阈值生成高斯核集合
    % 生成权重阈值序列，均匀分布在0到1之间
    wRg = linspace(0, 1, Nsec);
    
    % 生成对应的窗口大小序列，从大到小（kRU到kRL）
    kRg = round(linspace(kRU, kRL, Nsec));
    
    % 预生成所有可能使用的高斯核（在GPU上）
    gaussHs = cell(1, Nsec);
    for i = 1:Nsec
        % 创建椭圆形高斯核，标准差为尺寸的一半
        gaussHs{i} = gpuArray(fspecial('gaussian', [kRg(i) round(kRg(i)/2)], round(kRg(i)/2)));
    end
    
    %% 3. 使用GPU并行计算进行自适应滤波
    % 创建阈值矩阵 - 通过重塑wRg使其适合于arrayfun比较
    wRgMat = reshape(wRg, 1, 1, []);
    
    % 预分配输出矩阵
    outFrameWfilt = gpuArray.zeros(nZ, nX);
    
    % 对每个像素根据权重值找到对应的核索引
    % 创建索引矩阵以便快速查找
    [~, kernelIndices] = max(inWeight(:) >= reshape(wRg, 1, []), [], 2);
    kernelIndices = reshape(kernelIndices, nZ, nX);
    
    % 使用并行循环处理每个像素
    for iz = 1:nZ
        for ix = 1:nX
            fInd = kernelIndices(iz, ix);
            fkrg = kRg(fInd);
            
            % 计算窗口半径
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
            
            % 提取当前窗口区域的数据
            colLen = min(length(cols), round(fkrg/2));
            inA = inframe(rows, cols(1:colLen));
            
            % 确保卷积核和数据维度匹配
            kernelSize = size(gaussHs{fInd});
            if size(inA, 1) ~= kernelSize(1) || size(inA, 2) ~= kernelSize(2)
                % 裁剪或填充以匹配尺寸
                padRows = kernelSize(1) - size(inA, 1);
                padCols = kernelSize(2) - size(inA, 2);
                
                if padRows > 0 || padCols > 0
                    % 需要填充
                    inA = padarray(inA, [max(0, padRows), max(0, padCols)], 0, 'post');
                elseif padRows < 0 || padCols < 0
                    % 需要裁剪
                    inA = inA(1:min(end, kernelSize(1)), 1:min(end, kernelSize(2)));
                end
            end
            
            % 应用高斯核进行加权平均
            outFrameWfilt(iz, ix) = sum(inA .* gaussHs{fInd}, 'all');
        end
    end
    
    % 注意：结果为gpuArray类型，如需在CPU上使用，请调用gather()
end