%% PSOCT GPU加速处理脚本
% 此脚本展示如何使用GPU加速的PSOCT处理功能
% 作者：Yongxin Wang
% 日期：2023

% 添加函数路径
addpath('function');

% 设置数据路径和输出路径
data_path = ''; % 请设置您的数据路径
output_path = ''; % 请设置您的输出路径

% 检查数据路径
if isempty(data_path)
    data_path = uigetdir(pwd, '请选择数据文件夹');
    if data_path == 0
        error('没有选择数据文件夹，程序中止');
    end
end

% 检查输出路径
if isempty(output_path)
    output_path = uigetdir(pwd, '请选择输出文件夹');
    if output_path == 0
        error('没有选择输出文件夹，程序中止');
    end
end

% 创建输出文件夹（如果不存在）
if ~exist(output_path, 'dir')
    mkdir(output_path);
    fprintf('已创建输出文件夹: %s\n', output_path);
end

% 获取所有OCT文件
oct_files = dir(fullfile(data_path, '*.oct'));

if isempty(oct_files)
    error('在指定路径中找不到OCT文件: %s', data_path);
end

% 显示找到的文件
fprintf('找到 %d 个OCT文件:\n', length(oct_files));
for i = 1:length(oct_files)
    fprintf('[%d] %s\n', i, oct_files(i).name);
end

% 处理参数设置
do_avg = true;      % 是否执行平均处理
do_eig = false;     % 是否执行特征值处理
do_cfg1 = true;     % 是否使用配置1
do_cfg2 = false;    % 是否使用配置2
do_PhComp = true;   % 是否执行相位补偿
do_ssdopu = true;   % 是否使用split-spectrum DOPU
show_img = false;   % 是否显示图像
saveDicom = true;   % 是否保存为DICOM

% 处理文件
for i = 1:length(oct_files)
    % 获取完整文件路径
    file_path = fullfile(data_path, oct_files(i).name);
    [~, name, ~] = fileparts(oct_files(i).name);
    
    % 创建文件特定的输出目录
    file_output_path = fullfile(output_path, name);
    if ~exist(file_output_path, 'dir')
        mkdir(file_output_path);
    end
    
    % 显示处理信息
    fprintf('\n==================================================\n');
    fprintf('正在处理第 %d/%d 个文件: %s\n', i, length(oct_files), oct_files(i).name);
    
    % 调用GPU加速处理函数
    try
        rPSOCT_06_splitSpec_3_3D_gpu(...
            'file', file_path, ...
            'outDir', file_output_path, ...
            'display_name', name, ...
            'do_avg', do_avg, ...
            'do_eig', do_eig, ...
            'do_cfg1', do_cfg1, ...
            'do_cfg2', do_cfg2, ...
            'do_PhComp', do_PhComp, ...
            'do_ssdopu', do_ssdopu, ...
            'show_img', show_img, ...
            'saveDicom', saveDicom);
        
        fprintf('成功处理文件: %s\n', name);
    catch ME
        fprintf('处理文件 %s 时出错:\n%s\n', name, ME.message);
        fprintf('堆栈跟踪:\n');
        disp(getReport(ME, 'extended'));
    end
end

fprintf('\n==================================================\n');
fprintf('所有文件处理完成!\n');