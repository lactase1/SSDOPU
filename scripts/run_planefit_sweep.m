% run_planefit_sweep.m
% 该脚本用于扫描 Avnum 参数（从 3 到 30，步长 2），每次修改配置文件和主脚本中的相关参数，
% 然后运行主脚本，并记录每次运行的结果（包括时间、状态和输出目录）。
% 目的是评估不同 Avnum 值对平面拟合结果的影响。
%
% 使用方法：在 MATLAB 中运行此脚本（从工作区根目录或添加项目路径）。
% 注意：脚本会修改配置文件和主脚本，并在运行结束后恢复原始文件，以保持仓库清洁。

% 获取当前脚本文件的目录路径
root = fileparts(mfilename('fullpath'));
if isempty(root)
    root = pwd; % 如果从编辑器运行且无文件路径，则使用当前工作目录
end

% 构建项目根目录路径（脚本所在目录的父目录）
project_root = fullfile(root, '..');
project_root = matlab.lang.makeValidName(strrep(project_root, '\', '/'));

% 定义目标文件的绝对路径
cfg_file = fullfile('d:','\\1-Liu Jian','yongxin.wang','SSDOPU','function','config_params.m');
main_file = fullfile('d:','\\1-Liu Jian','yongxin.wang','SSDOPU','rPSOCT_06_splitSpec_3_3D_without_mat_wyx.m');

% 读取原始文件内容并备份
original_cfg = fileread(cfg_file);
original_main = fileread(main_file);

% 定义备份文件路径
backup_cfg = fullfile(fileparts(cfg_file),'config_params.m.bak');
backup_main = fullfile(fileparts(main_file),'rPSOCT_06_splitSpec_3_3D_without_mat_wyx.m.bak');

% 创建备份文件
fprintf('Backing up: %s -> %s\n', cfg_file, backup_cfg);
fid = fopen(backup_cfg, 'w'); fwrite(fid, original_cfg); fclose(fid);
fprintf('Backing up: %s -> %s\n', main_file, backup_main);
fid = fopen(backup_main, 'w'); fwrite(fid, original_main); fclose(fid);

% 初始化结果结构体和索引
results = struct(); % 用于存储每次运行的结果
idx = 0; % 结果索引

% 使用 try-catch 块处理全局错误，确保即使出错也能恢复文件
try
    % 循环扫描 Avnum 值：从 3 到 30，步长 2
    for Av = 3:2:30
        idx = idx + 1; % 递增结果索引

        % 为当前 Avnum 创建输出目录
        outdir = fullfile('D:','1-Liu Jian','yongxin.wang','tmp', sprintf('planefit_%dlayers', Av));
        if ~exist(outdir, 'dir')
            mkdir(outdir); % 如果目录不存在，则创建
        end

        % 打印当前运行信息
        fprintf('\n==== Run %d: Avnum=%d -> output: %s ====' , idx, Av, outdir);

        %% 步骤 1: 修改 config_params.m 文件中的 Avnum 参数
        % 读取配置文件内容
        cfg_text = fileread(cfg_file);
        % 使用正则表达式替换 Avnum 值
        cfg_text_new = regexprep(cfg_text, 'params\.polarization\.Avnum\s*=\s*\d+\s*;', sprintf('params.polarization.Avnum = %d;', Av));
        % 写入修改后的内容
        fid = fopen(cfg_file, 'w'); fwrite(fid, cfg_text_new); fclose(fid);
        fprintf('\n  updated %s (Avnum=%d)', cfg_file, Av);

        %% 步骤 2: 修改主脚本中的 output_base 行（基于行的安全替换，以保留反斜杠）
        % 读取主脚本内容
        main_text = fileread(main_file);
        % 将文本按行分割
        lines = regexp(main_text, '\n', 'split');
        found = false; % 标记是否找到 output_base 行
        for k = 1:length(lines)
            % 检查当前行是否匹配 output_base 赋值模式（忽略大小写和空格）
            if ~isempty(regexpi(strtrim(lines{k}), '^output_base\s*=\s*''.*''\s*;'))
                % 替换为新的输出目录
                lines{k} = sprintf('output_base = ''%s'';', outdir);
                found = true; % 标记已找到
                break; % 找到后退出循环
            end
        end
        if ~found
            % 如果未找到，则在文件顶部添加
            lines = [{sprintf("output_base = '%s';", outdir)}, lines];
        end
        % 将修改后的行重新组合成文本
        main_text_new = strjoin(lines, '\n');
        % 写入修改后的内容
        fid = fopen(main_file, 'w'); fwrite(fid, main_text_new); fclose(fid);
        fprintf('\n  updated %s output_base -> %s', main_file, outdir);

        %% 步骤 3: 运行主脚本并计时
        tstart = tic; % 开始计时
        try
            % 运行主脚本
            run(main_file);
            elapsed = toc(tstart); % 计算运行时间
            fprintf('\n  Run finished in %.1f s', elapsed);
            % 记录成功结果
            results(idx).Avnum = Av;
            results(idx).output = outdir;
            results(idx).time_s = elapsed;
            results(idx).status = 'ok';
        catch ME
            elapsed = toc(tstart); % 计算运行时间（即使失败）
            fprintf('\n  Run FAILED in %.1f s: %s', elapsed, ME.message);
            % 记录失败结果
            results(idx).Avnum = Av;
            results(idx).output = outdir;
            results(idx).time_s = elapsed;
            results(idx).status = 'error';
            results(idx).error = ME; % 保存错误信息
            % 继续下一个循环
        end

        % 运行间小暂停，避免系统过载
        pause(1);
    end
catch ME_global
    % 处理全局错误
    fprintf('\nSweep aborted due to error: %s\n', ME_global.message);
end

% 恢复原始文件，以保持仓库清洁
fprintf('\nRestoring original files...\n');
fid = fopen(cfg_file, 'w'); fwrite(fid, original_cfg); fclose(fid);
fid = fopen(main_file, 'w'); fwrite(fid, original_main); fclose(fid);

% 删除备份文件
delete(backup_cfg);
delete(backup_main);
fprintf('Backup files removed.\n');