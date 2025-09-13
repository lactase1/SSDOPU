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
%
% 示例:
%   % 将数据转移到GPU
%   gIMG_ch1 = gpuArray(IMG_ch1);
%   gIMG_ch2 = gpuArray(IMG_ch2);
%   % 在GPU上计算Stokes参数
%   [gS0, gS1, gS2, gS3] = cumulativeQUV_gpu(gIMG_ch1, gIMG_ch2);
%   % 如需，将结果传回CPU
%   S0 = gather(gS0);
%   S1 = gather(gS1);
%   S2 = gather(gS2);
%   S3 = gather(gS3);

    % 确保输入数据在GPU上
    if ~isa(IMG_ch1, 'gpuArray')
        IMG_ch1 = gpuArray(IMG_ch1);
    end
    
    if ~isa(IMG_ch2, 'gpuArray')
        IMG_ch2 = gpuArray(IMG_ch2);
    end
    
    % 计算相位差
    axis = angle(IMG_ch2 .* conj(IMG_ch1));
    
    % 在GPU上计算Stokes参数
    S0 = abs(IMG_ch1).^2 + abs(IMG_ch2).^2;
    S1 = abs(IMG_ch1).^2 - abs(IMG_ch2).^2;
    S2 = 2 .* abs(IMG_ch1) .* abs(IMG_ch2) .* cos(axis);
    S3 = 2 .* abs(IMG_ch1) .* abs(IMG_ch2) .* sin(-axis);
    
    % 注意：结果为gpuArray类型，如需在CPU上使用，请调用gather()
end