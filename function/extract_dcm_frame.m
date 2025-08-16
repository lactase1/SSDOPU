function extract_frame_from_dcm(input_folder, frame_number)
    % 从文件夹中的所有DCM文件中提取指定帧并保存为TIFF文件
    % 输入参数:
    %   input_folder: 包含DCM文件的文件夹路径
    %   frame_number: 要提取的帧号（从1开始计数），默认为301
    
    % 如果未提供参数，则使用默认值
    if nargin < 2
        frame_number = 301;
    end
    
    % 检查输入文件夹是否存在
    if ~exist(input_folder, 'dir')
        error(['指定的文件夹不存在: ' input_folder]);
    end
    
    % 获取所有DCM文件
    dcm_files = dir(fullfile(input_folder, '*.dcm'));
    
    if isempty(dcm_files)
        warning('指定路径中没有找到 .dcm 文件');
        return;
    end
    
    % 显示找到的文件数量
    fprintf('找到 %d 个 DCM 文件:\n', length(dcm_files));
    
    % 获取当前工作目录作为输出目录
    output_folder = pwd;
    
    % 处理每个DCM文件
    for i = 1:length(dcm_files)
        try
            % 构建完整文件路径
            dcm_file_path = fullfile(input_folder, dcm_files(i).name);
            fprintf('正在处理: %s\n', dcm_files(i).name);
            
            % 读取DCM文件
            ds = dicomread(dcm_file_path);
            
            % 检查是否为多帧数据
            if ndims(ds) < 3
                fprintf('警告: %s 不包含多帧数据，跳过\n', dcm_files(i).name);
                continue;
            end
            
            % 获取总帧数
            total_frames = size(ds, 3);
            fprintf('文件 %s 包含 %d 帧\n', dcm_files(i).name, total_frames);
            
            % 检查帧号是否有效
            if frame_number > total_frames
                fprintf('警告: %s 只有 %d 帧，无法提取第 %d 帧，跳过\n', ...
                    dcm_files(i).name, total_frames, frame_number);
                continue;
            end
            
            % 提取指定帧
            frame_data = ds(:, :, frame_number);
            
            % 生成输出文件名
            [~, base_name, ~] = fileparts(dcm_files(i).name);
            output_filename = sprintf('%s_frame%d.tiff', base_name, frame_number);
            output_path = fullfile(output_folder, output_filename);
            
            % 保存为TIFF文件
            imwrite(frame_data, output_path, 'Compression', 'none');
            fprintf('已保存: %s\n', output_filename);
            
        catch ME
            fprintf('处理 %s 时出错: %s\n', dcm_files(i).name, ME.message);
        end
    end
    
    fprintf('处理完成!\n');
    fprintf('输出文件保存在: %s\n', output_folder);
end

% 主函数 - 如果直接运行脚本则执行
if isequal(which('extract_frame_from_dcm'), mfilename('fullpath'))
    % 提示用户输入文件夹路径
    input_folder = input('请输入包含DCM文件的文件夹路径: ', 's');
    
    % 检查输入是否为空
    if isempty(input_folder)
        input_folder = pwd; % 使用当前目录作为默认值
        fprintf('使用当前目录: %s\n', input_folder);
    end
    
    % 调用处理函数
    extract_frame_from_dcm(input_folder, 301);
end