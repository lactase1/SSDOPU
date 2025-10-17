%% ========================================================================
% 脚本名: organize_results_by_oct.m
% 功能: 将不同参数配置的OCT结果按OCT数据重新组织
% 作者: GitHub Copilot
% 日期: 2025年9月24日
% 说明:
%   1. 遍历Output目录下的所有参数配置文件夹
%   2. 收集每个OCT数据的tiff文件
%   3. 以OCT数据名创建新文件夹，存储不同配置的结果
%   4. 支持选择特定类型的tiff文件进行组织
% ========================================================================

clc; clear; close all;

%% 配置参数
output_root = 'D:\1-Liu Jian\yongxin.wang\Output';
organized_root = 'D:\1-Liu Jian\yongxin.wang\Organized_Results';

% 选择要组织的tiff文件类型 (可以根据需要修改)
% 常见类型包括: 'Struc', 'dopu_SS', 'cumLA', 'drLA', 'PhR', 'Stokes' 等
selected_types = {
    '1-1_Struc_flat',           ... 结构图像(平坦化)
    '1-1_Struc',                ... 结构图像
    '1-3_1rep-Stokes',          ... 1次重复Stokes参数
    '1-3_4rep-Stokes',          ... 4次重复Stokes参数
    '1-4_dopu_SS',              ... DOPU图像
    '2-1_cumLA',                ... 累积线性退偏
    '2-2_cumLA_hsvColoring',    ... 累积线性退偏(HSV着色)
    '2-3_drLA',                 ... 动态范围线性退偏
    '2-4_drLA_hsvColoring',     ... 动态范围线性退偏(HSV着色)
    '2-5_PhR',                  ... 相位retardation
    '2-6_cumLA_rmBG',           ... 累积线性退偏(去背景)
    '2-7_cumLA_rmBG_hsvColoring', ... 累积线性退偏(去背景, HSV着色)
    '2-8_drLA_rmBG',            ... 动态范围线性退偏(去背景)
    '2-9_drLA_rmBG_hsvColoring', ... 动态范围线性退偏(去背景, HSV着色)
    '2-10_PhR_rmBG'             ... 相位retardation(去背景)
};

%% 预运行检查
fprintf('开始预运行检查...\n');

% 检查输出根目录
if ~exist(output_root, 'dir')
    error('输出根目录不存在: %s', output_root);
end

% 创建组织结果根目录
if ~exist(organized_root, 'dir')
    fprintf('创建组织结果根目录: %s\n', organized_root);
    mkdir(organized_root);
end

fprintf('预运行检查完成！\n\n');

%% 获取所有参数配置文件夹
fprintf('扫描参数配置文件夹...\n');
config_folders = dir(fullfile(output_root, 'ssdopu-*'));
config_folders = config_folders([config_folders.isdir]); % 只保留文件夹

fprintf('找到 %d 个参数配置文件夹:\n', length(config_folders));
for i = 1:length(config_folders)
    fprintf('  %d. %s\n', i, config_folders(i).name);
end
fprintf('\n');

%% 收集所有OCT数据
fprintf('收集所有OCT数据...\n');
oct_data_map = containers.Map(); % 用于存储OCT数据和对应配置的映射

for config_idx = 1:length(config_folders)
    config_name = config_folders(config_idx).name;
    config_path = fullfile(output_root, config_name);

    fprintf('处理配置 %d/%d: %s\n', config_idx, length(config_folders), config_name);

    % 获取该配置下的所有OCT文件夹
    oct_folders = dir(config_path);
    oct_folders = oct_folders([oct_folders.isdir] & ~strcmp({oct_folders.name}, '.') & ~strcmp({oct_folders.name}, '..'));

    for oct_idx = 1:length(oct_folders)
        oct_name = oct_folders(oct_idx).name;
        oct_path = fullfile(config_path, oct_name);
        tiff_path = fullfile(oct_path, 'tiff');

        % 检查tiff文件夹是否存在
        if ~exist(tiff_path, 'dir')
            fprintf('  警告: %s 的tiff文件夹不存在，跳过\n', oct_name);
            continue;
        end

        % 获取tiff文件
        tiff_files = dir(fullfile(tiff_path, '*.tiff'));

        % 初始化OCT数据映射
        if ~isKey(oct_data_map, oct_name)
            oct_data_map(oct_name) = struct();
        end

        % 存储该配置下的tiff文件
        oct_data_map(oct_name).(config_name) = tiff_files;

        fprintf('  %s: 找到 %d 个tiff文件\n', oct_name, length(tiff_files));
    end
end

fprintf('共收集到 %d 个OCT数据集\n\n', oct_data_map.length);

%% 组织结果
fprintf('开始组织结果...\n');
total_start_time = tic;

oct_names = keys(oct_data_map);
organized_count = 0;

for oct_idx = 1:length(oct_names)
    oct_name = oct_names{oct_idx};
    fprintf('组织OCT数据 %d/%d: %s\n', oct_idx, length(oct_names), oct_name);

    % 创建OCT结果文件夹
    oct_result_path = fullfile(organized_root, oct_name);
    if ~exist(oct_result_path, 'dir')
        mkdir(oct_result_path);
    end

    config_names = fieldnames(oct_data_map(oct_name));

    for config_idx = 1:length(config_names)
        config_name = config_names{config_idx};
        tiff_files = oct_data_map(oct_name).(config_name);

        % 创建配置子文件夹
        config_result_path = fullfile(oct_result_path, config_name);
        if ~exist(config_result_path, 'dir')
            mkdir(config_result_path);
        end

        % 筛选并复制选定的tiff文件
        copied_count = 0;
        for file_idx = 1:length(tiff_files)
            tiff_file = tiff_files(file_idx);
            file_name = tiff_file.name;

            % 检查是否为选定的类型
            is_selected = false;
            for type_idx = 1:length(selected_types)
                if contains(file_name, selected_types{type_idx})
                    is_selected = true;
                    break;
                end
            end

            if is_selected
                % 复制文件
                source_path = fullfile(tiff_files(file_idx).folder, file_name);
                dest_path = fullfile(config_result_path, file_name);

                try
                    copyfile(source_path, dest_path);
                    copied_count = copied_count + 1;
                catch ME
                    fprintf('    警告: 复制文件失败 %s -> %s: %s\n', source_path, dest_path, ME.message);
                end
            end
        end

        fprintf('  %s: 复制了 %d 个tiff文件\n', config_name, copied_count);
    end

    organized_count = organized_count + 1;
end

%% 完成总结
total_time = toc(total_start_time);
fprintf('\n=======================================\n');
fprintf('结果组织完成！\n');
fprintf('总耗时: %.1f秒 (%.1f分钟)\n', total_time, total_time/60);
fprintf('组织了 %d 个OCT数据集\n', organized_count);
fprintf('结果保存在: %s\n', organized_root);
fprintf('=======================================\n');

% 生成使用指南
fprintf('\n使用指南:\n');
fprintf('1. 每个OCT数据都有独立的文件夹\n');
fprintf('2. 每个OCT文件夹下包含不同参数配置的子文件夹\n');
fprintf('3. 每个配置文件夹包含选定的tiff结果文件\n');
fprintf('4. 支持的tiff类型: 结构图像、DOPU、Stokes参数、线性退偏、相位retardation等\n');
fprintf('5. 可以手动比较不同配置下相同OCT数据的效果\n');