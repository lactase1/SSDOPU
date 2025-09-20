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
data_path = 'D:\1-Liu Jian\yongxin.wang\PSOCT\tmp\';
output_root = 'D:\1-Liu Jian\yongxin.wang\PSOCT\2025-9-19\';

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
    % 第一组：从 4,9 开始的核心区域密集搜索 (1-30)
    4,  9;   % 组合1: 起始参数
    4, 11;   % 组合2
    4, 13;   % 组合3
    4, 15;   % 组合4
    4, 17;   % 组合5
    
    5,  7;   % 组合6
    5,  9;   % 组合7
    5, 11;   % 组合8: 之前最佳参数
    5, 13;   % 组合9
    5, 15;   % 组合10
    5, 17;   % 组合11
    
    6,  7;   % 组合12
    6,  9;   % 组合13
    6, 11;   % 组合14
    6, 13;   % 组合15
    6, 15;   % 组合16
    6, 17;   % 组合17
    
    7,  9;   % 组合18
    7, 11;   % 组合19
    7, 13;   % 组合20
    7, 15;   % 组合21
    7, 17;   % 组合22
    7, 19;   % 组合23
    
    8,  9;   % 组合24
    8, 11;   % 组合25
    8, 13;   % 组合26
    8, 15;   % 组合27
    8, 17;   % 组合28
    8, 19;   % 组合29
    8, 21;   % 组合30
    
    % 第二组：扩展到更大范围 (31-60)
    3,  5;   % 组合31
    3,  7;   % 组合32
    3,  9;   % 组合33
    3, 11;   % 组合34
    3, 13;   % 组合35
    3, 15;   % 组合36
    3, 17;   % 组合37
    
    9, 11;   % 组合38
    9, 13;   % 组合39
    9, 15;   % 组合40
    9, 17;   % 组合41
    9, 19;   % 组合42
    9, 21;   % 组合43
    9, 23;   % 组合44
    
    10, 13;  % 组合45
    10, 15;  % 组合46
    10, 17;  % 组合47
    10, 19;  % 组合48
    10, 21;  % 组合49
    10, 23;  % 组合50
    
    11, 15;  % 组合51
    11, 17;  % 组合52
    11, 19;  % 组合53
    11, 21;  % 组合54
    11, 23;  % 组合55
    11, 25;  % 组合56
    
    12, 15;  % 组合57
    12, 17;  % 组合58
    12, 19;  % 组合59
    12, 21;  % 组合60
    
    % 第三组：极值和边界探索 (61-100)
    1,  3;   % 组合61: 极细滤波
    1,  5;   % 组合62
    1,  7;   % 组合63
    1,  9;   % 组合64
    1, 11;   % 组合65
    1, 13;   % 组合66
    
    2,  3;   % 组合67
    2,  5;   % 组合68
    2,  7;   % 组合69
    2,  9;   % 组合70
    2, 11;   % 组合71
    2, 13;   % 组合72
    2, 15;   % 组合73
    
    13, 17;  % 组合74: 粗滤波
    13, 19;  % 组合75
    13, 21;  % 组合76
    13, 23;  % 组合77
    13, 25;  % 组合78
    
    14, 17;  % 组合79
    14, 19;  % 组合80
    14, 21;  % 组合81
    14, 23;  % 组合82
    14, 25;  % 组合83
    14, 27;  % 组合84
    
    15, 19;  % 组合85: 最粗滤波
    15, 21;  % 组合86
    15, 23;  % 组合87
    15, 25;  % 组合88
    15, 27;  % 组合89
    15, 29;  % 组合90
    
    % 补充测试组合 (91-100)
    4,  7;   % 组合91
    4, 19;   % 组合92
    5,  5;   % 组合93
    5, 19;   % 组合94
    6,  5;   % 组合95
    6, 19;   % 组合96
    7,  7;   % 组合97
    7, 21;   % 组合98
    8,  7;   % 组合99
    9,  9;   % 组合100
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
    
    % 构建新路径
    base_path = 'D:\1-Liu Jian\yongxin.wang\PSOCT\2025-9-19\';
    new_path = fullfile(base_path, param_name);
    new_path = strrep(new_path, '\', '\\'); % 转义反斜杠
    
    % 直接查找和替换第23行的output_base定义
    old_line = "output_base = 'D:\1-Liu Jian\yongxin.wang\PSOCT\2025-9-17\wyx';";
    new_line = sprintf("output_base = '%s';", new_path);
    
    % 执行替换
    new_content = strrep(content, old_line, new_line);
    
    % 验证替换是否成功
    if strcmp(content, new_content)
        % 如果内容没有改变，尝试更灵活的匹配
        fprintf('尝试更灵活的匹配...\n');
        
        % 使用正则表达式进行替换
        pattern = "output_base\s*=\s*'[^']*';";
        replacement = sprintf("output_base = '%s';", new_path);
        new_content = regexprep(content, pattern, replacement);
        
        if strcmp(content, new_content)
            error('输出路径替换失败，请检查主脚本中的output_base定义格式');
        end
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