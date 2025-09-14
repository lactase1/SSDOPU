function [LA_c_cfg2,PhR_c_cfg2,cumLA_cfg2] = calLAPhRcfg2_gpu(IMG_ch1,IMG_ch2,topLine,dopu_splitSpec,kRL,kRU,h1,h2,Avnum,wovWinF,gpu_available)
% CALLAPHRALL_GPU 使用GPU计算线性偏振特性和相位视网膜 (配置2)
%   [LA_c_cfg2,PhR_c_cfg2,cumLA_cfg2] = calLAPhRcfg2_gpu(IMG_ch1,IMG_ch2,topLine,dopu_splitSpec,kRL,kRU,h1,h2,Avnum,wovWinF,gpu_available)
%
%   此函数计算PSOCT图像的线性偏振特性和相位视网膜参数(配置2)，支持GPU加速
%
% 输入参数:
%   IMG_ch1       - 第一通道的复数数据
%   IMG_ch2       - 第二通道的复数数据
%   topLine       - 表面分割线，用于背景去除
%   dopu_splitSpec - 计算的dopu_splitSpec参数
%   kRL           - 高斯核最小尺寸
%   kRU           - 高斯核最大尺寸
%   h1            - 水平方向窗口大小1
%   h2            - 水平方向窗口大小2
%   Avnum         - 平均数量
%   wovWinF       - 是否使用无窗口滤波
%   gpu_available - GPU是否可用的标志
%
% 输出参数:
%   LA_c_cfg2     - 线性偏振特性 (配置2)
%   PhR_c_cfg2    - 相位视网膜 (配置2)
%   cumLA_cfg2    - 累积线性偏振特性 (配置2)

% 初始化输出变量
nZ = size(IMG_ch1, 1);
nX = size(IMG_ch1, 2);
LA_c_cfg2 = zeros(nZ-Avnum, nX, 3);
PhR_c_cfg2 = zeros(nZ-Avnum, nX);
cumLA_cfg2 = zeros(nZ-Avnum, nX, 3);

if gpu_available
    LA_c_cfg2 = gpuArray(LA_c_cfg2);
    PhR_c_cfg2 = gpuArray(PhR_c_cfg2);
    cumLA_cfg2 = gpuArray(cumLA_cfg2);
end

% 计算轴角
axis = angle(IMG_ch2.*conj(IMG_ch1));

% 计算Stokes参数
if gpu_available
    [S0,S1,S2,S3] = cumulativeQUV_gpu(IMG_ch1, IMG_ch2);
else
    % 内联实现cumulativeQUV功能
    S0=abs(IMG_ch1).^2+abs(IMG_ch2).^2;
    S1=abs(IMG_ch1).^2-abs(IMG_ch2).^2;
    S2=2.*abs(IMG_ch1).*abs(IMG_ch2).*cos(axis);
    S3=2.*abs(IMG_ch1).*abs(IMG_ch2).*sin(-axis);
end

% 计算规范化Stokes参数
Q = S1./S0;
U = S2./S0;
V = S3./S0;

% 计算局部线性偏振特性和相位视网膜 (配置2)
for i = 1:(nZ-Avnum)
    % 配置2 - 调整的滤波参数
    if wovWinF
        if gpu_available
            % 使用GPU加速的可变窗口高斯滤波 (配置2)
            tempAvgQ = vWinAvgFiltOpt_2_1_gpu(Q(i:i+Avnum-1,:), dopu_splitSpec(i:i+Avnum-1,:), kRL, kRU, h1);
            tempAvgU = vWinAvgFiltOpt_2_1_gpu(U(i:i+Avnum-1,:), dopu_splitSpec(i:i+Avnum-1,:), kRL, kRU, h1);
            tempAvgV = vWinAvgFiltOpt_2_1_gpu(V(i:i+Avnum-1,:), dopu_splitSpec(i:i+Avnum-1,:), kRL, kRU, h2);
        else
            % CPU版本的可变窗口高斯滤波 (配置2)
            tempAvgQ = vWinAvgFiltOpt_2_1(Q(i:i+Avnum-1,:), dopu_splitSpec(i:i+Avnum-1,:), kRL, kRU, h1);
            tempAvgU = vWinAvgFiltOpt_2_1(U(i:i+Avnum-1,:), dopu_splitSpec(i:i+Avnum-1,:), kRL, kRU, h1);
            tempAvgV = vWinAvgFiltOpt_2_1(V(i:i+Avnum-1,:), dopu_splitSpec(i:i+Avnum-1,:), kRL, kRU, h2);
        end
    else
        % 常规滤波
        tempAvgQ = mean(Q(i:i+Avnum-1,:), 1);
        tempAvgU = mean(U(i:i+Avnum-1,:), 1);
        tempAvgV = mean(V(i:i+Avnum-1,:), 1);
    end
    
    % 计算线性偏振特性和相位视网膜
    LA_c_cfg2(i,:,1) = tempAvgQ;
    LA_c_cfg2(i,:,2) = tempAvgU;
    LA_c_cfg2(i,:,3) = tempAvgV;
    PhR_c_cfg2(i,:) = sqrt(tempAvgQ.^2 + tempAvgU.^2);
    
    % 计算累积线性偏振特性
    if i == 1
        cumLA_cfg2(i,:,1) = tempAvgQ;
        cumLA_cfg2(i,:,2) = tempAvgU;
        cumLA_cfg2(i,:,3) = tempAvgV;
    else
        % 累积计算
        cumLA_cfg2(i,:,1) = cumLA_cfg2(i-1,:,1) + tempAvgQ;
        cumLA_cfg2(i,:,2) = cumLA_cfg2(i-1,:,2) + tempAvgU;
        cumLA_cfg2(i,:,3) = cumLA_cfg2(i-1,:,3) + tempAvgV;
    end
end

% 将结果转回CPU
if gpu_available
    LA_c_cfg2 = gather(LA_c_cfg2);
    PhR_c_cfg2 = gather(PhR_c_cfg2);
    cumLA_cfg2 = gather(cumLA_cfg2);
end

end