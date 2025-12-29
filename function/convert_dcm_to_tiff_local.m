%% ====================================================================================
% 函数名: convert_dcm_to_tiff_local
% 功能: 快速提取DCM文件指定帧并转换为单帧TIFF文件
% 输入参数:
%   dcm_folder - DCM文件所在文件夹路径
%   target_frame - 目标帧号(从1开始)
%   output_dir - 输出文件夹路径(可选，默认在 dcm_folder 下创建 tiff/png 子目录)
% 输出: 在指定目录内保存单帧TIFF和PNG文件
% ====================================================================================
function convert_dcm_to_tiff_local(dcm_folder, target_frame, output_dir)
    % 参数检查和默认值设置
    if nargin < 2
        error('至少需要提供DCM文件夹路径和目标帧号');
    end

    if nargin < 3 || isempty(output_dir)
        output_dir = dcm_folder;
    end

    % 检查DCM文件夹是否存在
    if ~exist(dcm_folder, 'dir')
        error(['DCM文件夹不存在: ' dcm_folder]);
    end

    % 输出目录：如果指定了 output_dir，则直接使用；否则在 dcm_folder 下创建子目录
    if strcmp(output_dir, dcm_folder)
        tiff_dir = fullfile(dcm_folder, 'tiff');
        png_dir  = fullfile(dcm_folder, 'png');
    else
        tiff_dir = output_dir;  % 直接使用传入的 tiff 目录
        png_dir  = fullfile(fileparts(output_dir), 'png');  % 在同级目录创建 png 文件夹
    end
    
    % 尝试创建子目录，如失败则回退到 DCM 文件夹本身
    try
        if ~exist(tiff_dir, 'dir'), mkdir(tiff_dir); end
    catch ME
        warning('无法创建 TIFF 输出目录 %s: %s\n将退回保存到 DCM 文件夹 %s', tiff_dir, ME.message, dcm_folder);
        tiff_dir = dcm_folder;
    end
    try
        if ~exist(png_dir, 'dir'), mkdir(png_dir); end
    catch ME
        warning('无法创建 PNG 输出目录 %s: %s\n将退回保存到 DCM 文件夹 %s', png_dir, ME.message, dcm_folder);
        png_dir = dcm_folder;
    end

    % 查找所有DCM文件
    dcm_files = dir(fullfile(dcm_folder, '*.dcm'));

    if isempty(dcm_files)
        warning(['在指定目录中未找到DCM文件: ' dcm_folder]);
        return;
    end

    fprintf('转换 %d 个DCM文件的第%d帧到 TIFF 和 PNG...\n', length(dcm_files), target_frame);

    % 处理每个DCM文件
    success_count = 0;
    for i = 1:length(dcm_files)
        dcm_filename = dcm_files(i).name;
        dcm_filepath = fullfile(dcm_folder, dcm_filename);
        
        try
            % 读取DCM文件
            dcm_data = dicomread(dcm_filepath);
            
            % 确定帧的位置
            if ndims(dcm_data) == 4
                [~, ~, ~, total_frames] = size(dcm_data);
                frame_dim = 4;
            elseif ndims(dcm_data) == 3
                [~, ~, total_frames] = size(dcm_data);
                frame_dim = 3;
            else
                % 2D图像，直接保存
                frame_dim = 0;
                total_frames = 1;
            end
            
            % 提取目标帧
            if frame_dim == 0
                frame_data = dcm_data;
            else
                actual_frame = min(target_frame, total_frames);
                if frame_dim == 4
                    frame_data = dcm_data(:, :, :, actual_frame);
                else
                    frame_data = dcm_data(:, :, actual_frame);
                end
            end
            
            % 生成输出文件名
            [~, base_name, ~] = fileparts(dcm_filename);
            tiff_filename = sprintf('%s_frame%d.tiff', base_name, target_frame);
            png_filename  = sprintf('%s_frame%d.png', base_name, target_frame);
            tiff_filepath = fullfile(tiff_dir, tiff_filename);
            png_filepath  = fullfile(png_dir, png_filename);

            % 尝试分别保存 TIFF 和 PNG
            saved_any = false;
            try
                imwrite(frame_data, tiff_filepath);
                fprintf('保存 TIFF: %s\n', tiff_filepath);
                saved_any = true;
            catch ME
                fprintf('保存 TIFF 失败 %s: %s\n', tiff_filepath, ME.message);
            end
            try
                imwrite(frame_data, png_filepath);
                fprintf('保存 PNG: %s\n', png_filepath);
                saved_any = true;
            catch ME
                fprintf('保存 PNG 失败 %s: %s\n', png_filepath, ME.message);
            end

            if saved_any
                success_count = success_count + 1;
            end
            
        catch ME
            fprintf('处理文件 %s 出错: %s\n', dcm_filename, ME.message);
            continue;
        end
    end

    fprintf('完成! 成功转换 %d/%d 个文件。输出目录: %s 以及 %s\n', success_count, length(dcm_files), tiff_dir, png_dir);
end
