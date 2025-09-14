function [S0,S1,S2,S3] = cumulativeQUV_gpu(IMG_ch1, IMG_ch2)
%% CUMULATIVEQUV_GPU 使用GPU计算Stokes参数
%   [S0,S1,S2,S3] = cumulativeQUV_gpu(IMG_ch1, IMG_ch2)
%
%   此函数计算基于两个通道复数数据的Stokes参数
%   使用GPU加速计算过程，显著提高大数据集处理速度
%
% 输入参数:
%   IMG_ch1 - 第一个通道的复数数据，可以是gpuArray
%   IMG_ch2 - 第二个通道的复数数据，可以是gpuArray
%
% 输出参数:
%   S0, S1, S2, S3 - 计算得到的Stokes参数，为gpuArray类型

% 计算相位差
axis = angle(IMG_ch2.*conj(IMG_ch1));

% 计算Stokes参数
S0 = abs(IMG_ch1).^2 + abs(IMG_ch2).^2;
S1 = abs(IMG_ch1).^2 - abs(IMG_ch2).^2;
S2 = 2.*abs(IMG_ch1).*abs(IMG_ch2).*cos(axis);
S3 = 2.*abs(IMG_ch1).*abs(IMG_ch2).*sin(-axis);

end