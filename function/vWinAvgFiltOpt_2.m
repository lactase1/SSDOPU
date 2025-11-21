function [outFrameWfilt] = vWinAvgFiltOpt_2(inFrame,inWeight,kRL,kRU,Nsec)
    % vWinAvgFiltOpt_2: 动态加权高斯滤波函数
    % 该函数根据输入权重动态选择高斯核大小，对输入数据进行加权滤波。
    % 适用于需要根据置信度调整滤波强度的场景，如偏振光学数据处理。
    %
    % 输入参数:
    %   inFrame - 输入数据矩阵 [nZ x nX]，待滤波的数据
    %   inWeight - 权重矩阵 [nZ x nX]，每个像素的置信度 (0-1)
    %   kRL - 高斯核大小的下限 (整数)
    %   kRU - 高斯核大小的上限 (整数)
    %   Nsec - 分段数 (可选，默认 kRU - kRL + 1)，控制高斯核变化的精细度
    %
    % 输出参数:
    %   outFrameWfilt - 滤波后的输出数据矩阵 [nZ x nX]

    % 处理可选参数 Nsec
    if nargin == 4, Nsec = kRU - kRL +1; end

    % 获取输入矩阵的尺寸
    [nZ,nX] = size(inFrame);

    % 计算扩展边界的大小 (基于最大核大小)
    krs = (kRU-1)/2; kcs = (kRU-1)/2;

    % 如果核大小为0，直接返回加权结果
    if krs == 0 && kcs ==0, outFrameWfilt = inFrame.*inWeight; return; end

    % 扩展输入数据和权重矩阵 (padding zeros)
    inframe = zeros(nZ+krs*2,nX+kcs*2); inweight = inframe;
    inframe(krs+1:end-krs,kcs+1:end-kcs) = inFrame;
    inweight(krs+1:end-krs,kcs+1:end-kcs) = inWeight;

    % 生成权重范围和对应的核大小
    wRg = linspace(0,1,Nsec);  % 权重范围 [0,1]
    kRg = round(linspace(kRU,kRL,Nsec));  % 核大小范围 [kRU, kRL]

    % 生成一系列高斯核
    for i = 1:1:Nsec
        gaussHs{i} = fspecial('gaussian',[kRg(i) kRg(i)],round(kRg(i)/2));
    end

    % 遍历每个像素，进行动态加权滤波
    for iz = 1:nZ
        for ix = 1:nX
            % 根据当前像素的权重选择合适的核索引
            wInds = wRg >= inWeight(iz,ix);
            pInds = find(wInds);
            fInd = pInds(1);  % 选择第一个满足条件的索引
            fkrg = kRg(fInd);  % 对应的核大小

            % 计算核的半径范围
            if mod(fkrg,2)
                ks1 = (fkrg-1)/2; ks2 = ks1;
            else
                ks1 = fkrg/2-1;
                ks2 = ks1+1;
            end

            % 定义局部区域的行列索引
            cols = ix+kcs-ks1:ix+kcs+ks2;
            rows = iz+krs-ks1:iz+krs+ks2;

            % 提取局部数据
            inA = inframe(rows,cols);

            % 使用选定的高斯核进行加权滤波
            outFrameWfilt(iz,ix) = sum(inA.*gaussHs{fInd},'all');
        end
    end
end