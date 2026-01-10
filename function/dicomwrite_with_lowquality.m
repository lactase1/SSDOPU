function dicomwrite_with_lowquality(data, dcm_path, low_quality_dir, scale_factor)
% DICOMWRITE_WITH_LOWQUALITY 同时保存原始DCM和低质量压缩版本
%
% 输入参数:
%   data            - 要保存的图像数据 (uint8)
%   dcm_path        - 原始DCM文件的完整路径
%   low_quality_dir - 低质量DCM文件保存的目录（空字符串则跳过低质量版本生成）
%   scale_factor    - 缩放因子 (0-1之间，默认0.5表示缩小到50%)
%
% 功能:
%   1. 保存原始DCM文件
%   2. 对图像进行下采样后保存到low_quality_dir（如果启用）

    if nargin < 4
        scale_factor = 0.5;  % 默认缩小到50%
    end
    
    % 1. 保存原始DCM文件
    dicomwrite(data, dcm_path);
    
    % 2. 如果low_quality_dir为空或scale_factor<=0，则跳过低质量版本生成
    if isempty(low_quality_dir) || scale_factor <= 0
        return;
    end
    
    % 3. 创建低质量目录（如果不存在）
    if ~exist(low_quality_dir, 'dir')
        mkdir(low_quality_dir);
    end
    
    % 4. 生成低质量版本
    try
        % 获取文件名
        [~, filename, ext] = fileparts(dcm_path);
        low_quality_path = fullfile(low_quality_dir, [filename, ext]);
        
        % 获取数据维度
        dims = size(data);
        ndims_data = length(dims);
        
        if ndims_data == 2
            % 2D 灰度图像 [H, W]
            new_size = round([dims(1), dims(2)] * scale_factor);
            data_low = imresize(data, new_size, 'bilinear');
            
        elseif ndims_data == 3
            % 3D 数据：可能是 [H, W, C] RGB 或 [H, W, Frames] 多帧灰度
            if dims(3) == 3
                % RGB 图像 [H, W, 3]
                new_size = round([dims(1), dims(2)] * scale_factor);
                data_low = imresize(data, new_size, 'bilinear');
            else
                % 多帧灰度 [H, W, Frames]
                new_h = max(1, round(dims(1) * scale_factor));
                new_w = max(1, round(dims(2) * scale_factor));
                data_low = zeros(new_h, new_w, dims(3), class(data));
                for f = 1:dims(3)
                    data_low(:, :, f) = imresize(data(:, :, f), [new_h, new_w], 'bilinear');
                end
            end
            
        elseif ndims_data == 4
            % 4D 数据：[H, W, C, Frames] 或 [H, W, 1, Frames]
            new_h = max(1, round(dims(1) * scale_factor));
            new_w = max(1, round(dims(2) * scale_factor));
            data_low = zeros(new_h, new_w, dims(3), dims(4), class(data));
            
            for f = 1:dims(4)
                if dims(3) == 1
                    % 灰度多帧 [H, W, 1, Frames]
                    data_low(:, :, 1, f) = imresize(squeeze(data(:, :, 1, f)), [new_h, new_w], 'bilinear');
                elseif dims(3) == 3
                    % RGB 多帧 [H, W, 3, Frames]
                    for c = 1:3
                        data_low(:, :, c, f) = imresize(squeeze(data(:, :, c, f)), [new_h, new_w], 'bilinear');
                    end
                else
                    % 其他情况，逐通道处理
                    for c = 1:dims(3)
                        data_low(:, :, c, f) = imresize(squeeze(data(:, :, c, f)), [new_h, new_w], 'bilinear');
                    end
                end
            end
        else
            % 不支持的维度，直接复制原始数据
            fprintf('警告: 不支持的数据维度 (%d)，低质量版本将与原始相同\\n', ndims_data);
            data_low = data;
        end
        
        % 保存低质量DCM
        dicomwrite(data_low, low_quality_path);
        
    catch ME
        fprintf('生成低质量DCM时出错: %s\\n', ME.message);
        % 出错时不影响主流程
    end
end
