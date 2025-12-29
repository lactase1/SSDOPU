function overwrite_config_params(updates, config_file)
% overwrite_config_params  简单直接地替换 config_params.m 中的指定赋值行（不备份）
%   updates: 结构体，支持字段示例：
%     updates.polarization.Avnum
%     updates.filters.h2_sigma
%     updates.filters.h2_size  (% [sz sz])
%     updates.mode.wovWinF
%   config_file: config_params.m 的完整路径
%
% 注意：此函数不会创建备份，也不会做复杂的验证，仅按正则替换文本并写回文件。

if nargin < 2
    error('overwrite_config_params: 需要两个参数: updates, config_file');
end
if ~isstruct(updates)
    error('overwrite_config_params: updates 必须为结构体');
end
if ~exist(config_file, 'file')
    error('overwrite_config_params: 找不到配置文件: %s', config_file);
end

% 读取文件内容
txt = fileread(config_file);

% 记录替换结果
repl_report = struct();

% 替换 Avnum
if isfield(updates, 'polarization') && isfield(updates.polarization, 'Avnum')
    pattern = 'params\.polarization\.Avnum\s*=\s*\d+\s*;';
    matches = regexp(txt, pattern, 'match');
    repl_report.Avnum_matches_before = numel(matches);
    if repl_report.Avnum_matches_before > 0
        txt = regexprep(txt, pattern, sprintf('params.polarization.Avnum = %d;', updates.polarization.Avnum));
        repl_report.Avnum_replaced = true;
    else
        repl_report.Avnum_replaced = false;
    end
end

% 替换 h2_sigma
if isfield(updates, 'filters') && isfield(updates.filters, 'h2_sigma')
    pattern = 'params\.filters\.h2_sigma\s*=\s*[\d\.Ee\+\-]+\s*;';
    matches = regexp(txt, pattern, 'match');
    repl_report.h2_sigma_matches_before = numel(matches);
    if repl_report.h2_sigma_matches_before > 0
        txt = regexprep(txt, pattern, sprintf('params.filters.h2_sigma = %g;', updates.filters.h2_sigma));
        repl_report.h2_sigma_replaced = true;
    else
        repl_report.h2_sigma_replaced = false;
    end
end

% 替换 h2_size
if isfield(updates, 'filters') && isfield(updates.filters, 'h2_size')
    sz = updates.filters.h2_size;
    pattern = 'params\.filters\.h2_size\s*=\s*\[[^\]]*\]\s*;';
    matches = regexp(txt, pattern, 'match');
    repl_report.h2_size_matches_before = numel(matches);
    if repl_report.h2_size_matches_before > 0
        if numel(sz) == 2
            txt = regexprep(txt, pattern, sprintf('params.filters.h2_size = [%d %d];', round(sz(1)), round(sz(2))));
        else
            txt = regexprep(txt, pattern, sprintf('params.filters.h2_size = [%s];', num2str(sz)));
        end
        repl_report.h2_size_replaced = true;
    else
        repl_report.h2_size_replaced = false;
    end
end

% 替换 wovWinF
if isfield(updates, 'mode') && isfield(updates.mode, 'wovWinF')
    pattern = 'params\.mode\.wovWinF\s*=\s*\d+\s*;';
    matches = regexp(txt, pattern, 'match');
    repl_report.wovWinF_matches_before = numel(matches);
    if repl_report.wovWinF_matches_before > 0
        txt = regexprep(txt, pattern, sprintf('params.mode.wovWinF = %d;', updates.mode.wovWinF));
        repl_report.wovWinF_replaced = true;
    else
        repl_report.wovWinF_replaced = false;
    end
end

% 写回文件（覆盖）
fid = fopen(config_file, 'w');
if fid == -1
    error('overwrite_config_params: 无法打开文件写入: %s', config_file);
end
fwrite(fid, txt);
fclose(fid);

% 读取并打印替换后的证据（行片段）
new_txt = fileread(config_file);
if isfield(repl_report, 'Avnum_replaced')
    if repl_report.Avnum_replaced
        fprintf('overwrite_config_params: Avnum 替换成功，原匹配数=%d\n', repl_report.Avnum_matches_before);
        m = regexp(new_txt, 'params\.polarization\.Avnum\s*=\s*\d+\s*;', 'match');
        fprintf('新行: %s\n', m{1});
    else
        warning('overwrite_config_params: 未找到 Avnum 赋值行以进行替换');
    end
end
if isfield(repl_report, 'h2_sigma_replaced')
    if repl_report.h2_sigma_replaced
        fprintf('overwrite_config_params: h2_sigma 替换成功，原匹配数=%d\n', repl_report.h2_sigma_matches_before);
        m = regexp(new_txt, 'params\.filters\.h2_sigma\s*=\s*[\d\.Ee\+\-]+\s*;', 'match');
        fprintf('新行: %s\n', m{1});
    else
        warning('overwrite_config_params: 未找到 h2_sigma 赋值行以进行替换');
    end
end
if isfield(repl_report, 'h2_size_replaced')
    if repl_report.h2_size_replaced
        fprintf('overwrite_config_params: h2_size 替换成功，原匹配数=%d\n', repl_report.h2_size_matches_before);
        m = regexp(new_txt, 'params\.filters\.h2_size\s*=\s*\[[^\]]*\]\s*;', 'match');
        fprintf('新行: %s\n', m{1});
    else
        warning('overwrite_config_params: 未找到 h2_size 赋值行以进行替换');
    end
end
if isfield(repl_report, 'wovWinF_replaced')
    if repl_report.wovWinF_replaced
        fprintf('overwrite_config_params: wovWinF 替换成功，原匹配数=%d\n', repl_report.wovWinF_matches_before);
        m = regexp(new_txt, 'params\.mode\.wovWinF\s*=\s*\d+\s*;', 'match');
        fprintf('新行: %s\n', m{1});
    else
        warning('overwrite_config_params: 未找到 wovWinF 赋值行以进行替换');
    end
end

fprintf('已直接写回配置到 %s\n', config_file);
end