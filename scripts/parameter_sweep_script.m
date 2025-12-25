%% ========================================================================
% 脚本名: parameter_sweep_script.m
% 功能: 自动修改参数并运行PS-OCT处理脚本，测试不同参数组合的效果
% 作者: GitHub Copilot
% 日期: 2025年9月19日
% 说明: 
%   1. 直接修改 config_params.m 文件中的参数值（最简单有效的方法）
%   2. 修改主处理脚本中的输出路径为 ssdopu+参数
%   3. 运行主处理脚本 rPSOCT_06_splitSpec_3_3D_without_mat_wyx.m
%   4. 重复以上步骤100次，测试100组不同参数
% ========================================================================

clc; clear; close all;

%% 文件路径设置
config_file = 'D:\1-Liu Jian\yongxin.wang\SSDOPU\function\config_params.m';
main_script = 'D:\1-Liu Jian\yongxin.wang\SSDOPU\rPSOCT_06_splitSpec_3_3D_without_mat_wyx.m';
data_path = '/d/1-Liu Jian/yongxin.wang/nail';
output_root = 'D:\1-Liu Jian\yongxin.wang\Output';

%% 预运行检查
fprintf('开始预运行检查...\n');

% 检查文件是否存在
if ~exist(config_file, 'file')
    error('配置文件不存在: %s', config_file);
end
if ~exist(main_script, 'file')
    error('主处理脚本不存在: %s', main_script);
end

% 检查数据路径
if ~exist(data_path, 'dir')
    warning('数据路径不存在: %s', data_path);
    fprintf('尝试创建数据目录...\n');
    try
        mkdir(data_path);
        fprintf('已创建数据目录: %s\n', data_path);
        warning('请确保将数据文件放入此目录后再运行脚本');
    catch
        error('无法创建数据目录: %s', data_path);
    end
else
    fprintf('数据路径检查通过: %s\n', data_path);
end

% 检查和创建输出根目录
if ~exist(output_root, 'dir')
    fprintf('输出根目录不存在，正在创建: %s\n', output_root);
    try
        mkdir(output_root);
        fprintf('已创建输出根目录: %s\n', output_root);
    catch
        error('无法创建输出根目录: %s', output_root);
    end
else
    fprintf('输出根目录检查通过: %s\n', output_root);
end

fprintf('预运行检查完成！\n\n');

%% 参数设置
% 定义100组不同的滤波核参数组合
% 重要约束：kRU_cfg1 必须为奇数！
% 从 4,9 开始，围绕 5,11 进行全面参数扫描
param_combinations = [
    %% 刘老师要求参数
    % 3 21;
    % 5 21;
    % 7 21;
    9 21;
];

fprintf('开始参数扫描...\n');
fprintf('将测试 %d 组参数组合\n', size(param_combinations, 1));

%% 备份原始文件
backup_config = [config_file '.backup'];
backup_main = [main_script '.backup'];

copyfile(config_file, backup_config);
copyfile(main_script, backup_main);
fprintf('已备份原始文件\n');

%% 参数扫描主循环
total_start_time = tic;

for param_idx = 1:size(param_combinations, 1)
    % 获取当前参数组合
    kRL = param_combinations(param_idx, 1);
    kRU = param_combinations(param_idx, 2);
    param_name = sprintf('ssdopu-kRL_%02d-kRU_%02d', kRL, kRU);
    
    fprintf('\n========================================\n');
    fprintf('开始测试参数组合 %d/%d: %s\n', param_idx, size(param_combinations, 1), param_name);
    fprintf('kRL_cfg1=%d, kRU_cfg1=%d\n', kRL, kRU);
    fprintf('========================================\n');
    
    try
        % 步骤0: 创建当前参数组合的输出目录
        current_output_dir = fullfile(output_root, param_name);
        if ~exist(current_output_dir, 'dir')
            mkdir(current_output_dir);
            fprintf('0. 创建输出目录: %s\n', current_output_dir);
        end
        
        % 步骤1: 修改配置文件中的参数
        fprintf('1. 修改配置文件参数...\n');
        modify_config_params(config_file, kRL, kRU);
        
        % 步骤2: 修改主脚本中的输出路径
        fprintf('2. 修改输出路径...\n');
        modify_output_path(main_script, param_name);
        
        % 步骤3: 运行主处理脚本
        fprintf('3. 运行处理脚本...\n');
        param_start_time = tic;
        
        % 运行脚本
        run(main_script);
        
        param_time = toc(param_start_time);
        fprintf('参数组合 %d 完成，耗时: %.1f秒 (%.1f分钟)\n', param_idx, param_time, param_time/60);
        
    catch ME
        fprintf('处理参数组合 %d 时出错: %s\n', param_idx, ME.message);
        if ~isempty(ME.stack)
            fprintf('错误位置: %s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
        end
        continue; % 继续处理下一组参数
    end
end

%% 恢复原始文件
fprintf('\n恢复原始文件...\n');
copyfile(backup_config, config_file);
copyfile(backup_main, main_script);
delete(backup_config);
delete(backup_main);

%% 扫描完成总结
total_time = toc(total_start_time);
fprintf('\n=======================================\n');
fprintf('参数扫描完成！\n');
fprintf('总耗时: %.1f秒 (%.1f分钟)\n', total_time, total_time/60);
fprintf('测试了 %d 组参数组合\n', size(param_combinations, 1));
fprintf('=======================================\n');

% 生成结果比较指南
fprintf('\n结果比较建议:\n');
fprintf('1. 检查各个 ssdopu-kRL_xx-kRU_xx 文件夹中的处理结果\n');
fprintf('2. 重点查看 DOPU 图像 (_1-4_dopu_SS.dcm) 的质量\n');
fprintf('3. 比较结构图像 (_1-1_Struc.dcm) 的清晰度\n');
fprintf('4. 选择视觉效果最佳的参数组合\n');

%% ========================================================================
% 函数名: modify_config_params
% 功能: 直接修改 config_params.m 文件中的参数值
% ========================================================================
function modify_config_params(config_file, kRL_cfg1, kRU_cfg1)
    % 读取文件内容
    content = fileread(config_file);
    
    % 使用正则表达式替换参数值
    % 替换 kRL_cfg1 参数
    pattern1 = 'params\.polarization\.kRL_cfg1\s*=\s*\d+\s*;';
    replacement1 = sprintf('params.polarization.kRL_cfg1 = %d;', kRL_cfg1);
    content = regexprep(content, pattern1, replacement1);
    
    % 替换 kRU_cfg1 参数
    pattern2 = 'params\.polarization\.kRU_cfg1\s*=\s*\d+\s*;';
    replacement2 = sprintf('params.polarization.kRU_cfg1 = %d;', kRU_cfg1);
    content = regexprep(content, pattern2, replacement2);
    
    % 写回文件
    fid = fopen(config_file, 'w');
    if fid == -1
        error('无法写入配置文件: %s', config_file);
    end
    fwrite(fid, content, 'char');
    fclose(fid);
    
    % 验证修改是否成功
    validate_config_modification(config_file, kRL_cfg1, kRU_cfg1);
    
    fprintf('  更新配置文件: kRL_cfg1=%d, kRU_cfg1=%d\n', kRL_cfg1, kRU_cfg1);
end

%% ========================================================================
% 函数名: validate_config_modification
% 功能: 验证配置文件参数修改是否成功
% ========================================================================
function validate_config_modification(config_file, expected_kRL, expected_kRU)
    try
        % 读取修改后的文件内容
        content = fileread(config_file);
        
        % 提取实际的参数值
        kRL_match = regexp(content, 'params\.polarization\.kRL_cfg1\s*=\s*(\d+)', 'tokens');
        kRU_match = regexp(content, 'params\.polarization\.kRU_cfg1\s*=\s*(\d+)', 'tokens');
        
        if ~isempty(kRL_match) && ~isempty(kRU_match)
            actual_kRL = str2double(kRL_match{1}{1});
            actual_kRU = str2double(kRU_match{1}{1});
            
            if actual_kRL ~= expected_kRL || actual_kRU ~= expected_kRU
                error('参数修改验证失败: 预期 kRL=%d, kRU=%d, 实际 kRL=%d, kRU=%d', ...
                    expected_kRL, expected_kRU, actual_kRL, actual_kRU);
            end
        else
            error('无法在配置文件中找到参数');
        end
    catch ME
        error('参数修改验证失败: %s', ME.message);
    end
end

%% ========================================================================
% 函数名: modify_output_path
% 功能: 修改主脚本中的输出路径
% ========================================================================
function modify_output_path(main_script, param_name)
    % 读取文件内容
    content = fileread(main_script);
    
    % 标准化路径处理 - 移除末尾的反斜杠以避免fullfile重复
    base_path = 'D:\1-Liu Jian\yongxin.wang\Output';
    new_path = fullfile(base_path, param_name);
    
    % 直接查找和替换当前已知的output_base行
    old_line = "output_base = 'D:\1-Liu Jian\yongxin.wang\Output';";
    new_line = sprintf("output_base = '%s';", new_path);
    
    % 执行替换
    new_content = strrep(content, old_line, new_line);
    
    % 如果没有找到匹配，尝试其他可能的格式
    if strcmp(content, new_content)
        fprintf('  尝试其他格式...\n');
        
        % 尝试没有反斜杠的格式（如果文件中有错误格式）
        old_line2 = "output_base = 'D:1-Liu Jianyongxin.wangOutput';";
        new_content = strrep(content, old_line2, new_line);
        
        if strcmp(content, new_content)
            error('无法找到可替换的output_base行，请手动检查主脚本');
        end
    end
    
    % 验证修改是否成功
    if ~contains(new_content, new_path)
        error('路径替换验证失败：新内容不包含期望路径 %s', new_path);
    end
    
    % 写回文件
    fid = fopen(main_script, 'w');
    if fid == -1
        error('无法写入主脚本: %s', main_script);
    end
    fwrite(fid, new_content, 'char');
    fclose(fid);
    
    fprintf('  更新输出路径: %s\n', new_path);
end