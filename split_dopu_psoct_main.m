%%% rPSOCT_process_single_file.m
%%% 1. frame based code to cal 3D LAs and Phase[two cfg: cfg1== LA, cfg2==Phase]
%%% 1. use doup to apply a variable gaussian fitler to Q U V before cal LA
%%%   and phase
%%% 2. speed the cal process, params.range.setZrg = 100; params.range.setZrg = 0 for process all
%%% ______________________________20240913_______________________________________________
%%% 3. optional to toggle vWinF
%%% 4.cfg2 was optinal to [do SVD with remove bias] and DDG
%%% ______________________________20240721_______________________________________________
%%% 5. optinal to perform split-spectrum DOPU
%%% 6. DIR: 03_2024.07.03.psOCTry2


% 添加function文件夹到搜索路径
script_dir = fileparts(mfilename('fullpath'));
function_path = fullfile(script_dir, 'function');
if exist(function_path, 'dir')
    addpath(function_path);
end

% 设置数据路径
data_path   = 'C:\yongxin.wang\glaucoma';
% σ * 6 + 1 // σ * 4 + 1
output_base = 'D:\1-Liu Jian\yongxin.wang\Output\G\ddg_3layer_3_3';
if ~exist(data_path, 'dir')
    error(['数据路径不存在: ' data_path]);
end

% 创建输出路径
if ~exist(output_base, 'dir')
    fprintf(['输出路径不存在: ' output_base]);
    mkdir(output_base);
end

% 获取所有oct文件
disp('正在处理数据...');
oct_files = dir(fullfile(data_path, '*.oct'));

if isempty(oct_files)
else
    % 显示找到的文件
    fprintf('找到 %d 个 OCT 文件:\n', length(oct_files));
    for i = 1:length(oct_files)
        fprintf('[%d] %s\n', i, oct_files(i).name);
    end
    
    % 在批处理模式下创建本次运行的根输出目录（便于管理多个文件的输出）
    if length(oct_files) > 1
        run_tag = datestr(now,'yyyymmdd_HHMMSS');
        run_output_base = fullfile(output_base, ['run_' run_tag]);
        if ~exist(run_output_base, 'dir')
            mkdir(run_output_base);
        end
        fprintf('创建本次批处理输出根目录: %s\n', run_output_base);
    else
        run_output_base = output_base;
    end

    % 逐个处理文件
    % 如果是批处理模式，开始整体计时以显示总进度/ETA
    total_files = length(oct_files);
    runTimer = tic;
    % 初始化总进度（每文件细分为100个单位，便于显示文件内的进度细分）
    print_progress(0, total_files*100, '总进度', 40, 'blocks', NaN);
    for i = 1:total_files
        full_path = fullfile(data_path, oct_files(i).name);
        fprintf('\n==================================================\n');
        fprintf('开始处理第 %d/%d 个文件: %s\n', i, total_files, oct_files(i).name);
        fprintf('==================================================\n');
        
        try
            % 记录本文件开始时间
            file_timer = tic;
            
            % 调用单文件处理函数（传入本次运行的根输出目录 run_output_base）
            % 传递文件索引、总文件数与总运行计时器以便内部更新总进度条
            rPSOCT_process_single_file(full_path, run_output_base, i, total_files, runTimer);
            
            % 计算并显示处理时间
            proc_time = toc(file_timer);
            fprintf('文件 %s 处理成功, 耗时: %.2f 秒 (%.2f 分钟)\n', oct_files(i).name, proc_time, proc_time/60);
            
            % 更新总进度（以每文件为 100 单位的细分单位），并显示 ETA
            globalTotal = total_files * 100;
            globalCurrent = i * 100; % 完成当前文件
            elapsed = toc(runTimer);
            if globalCurrent > 0
                avg_unit = elapsed / max(1, globalCurrent);
                eta = avg_unit * max(0, globalTotal - globalCurrent);
            else
                eta = NaN;
            end
            prefix_total = sprintf('总进度 %d/%d', i, total_files);
            print_progress(globalCurrent, globalTotal, prefix_total, 40, 'blocks', eta);
        catch ME
            fprintf('处理文件 %s 时出错: %s\n', oct_files(i).name, ME.message);
            % 将该文件计为已完成以推进总进度（避免停滞）
            globalTotal = total_files * 100;
            globalCurrent = i * 100;
            elapsed = toc(runTimer);
            if globalCurrent > 0
                avg_unit = elapsed / max(1, globalCurrent);
                eta = avg_unit * max(0, globalTotal - globalCurrent);
            else
                eta = NaN;
            end
            prefix_total = sprintf('总进度 %d/%d (含失败)', i, total_files);
            print_progress(globalCurrent, globalTotal, prefix_total, 40, 'blocks', eta);
            % 继续处理下一个文件
        end
    end
    
    fprintf('\n所有文件处理完成!\n');
    
end


function rPSOCT_process_single_file(varargin)
    clc;
    %% (1) 获取输入参数
    if nargin < 1
        error('请提供输入文件路径作为参数');
    end

    input_file_path = varargin{1};

    % 如果提供了输出路径，则使用提供的路径，否则使用默认路径
    if nargin >= 2
        output_base = varargin{2};
    else
        output_base = 'D:1-Liu Jianyongxin.wang	mpplanefit_21layers';
    end

    % 检查输入文件是否存在
    if ~exist(input_file_path, 'file')
        error(['输入文件不存在: ' input_file_path]);
    end

    % 创建输出目录
    if ~exist(output_base, 'dir')
        mkdir(output_base) 
    end

    % 可选的整体进度追踪参数（fileIdx, totalFiles, runTimer）
    fileIdx = [];
    totalFiles = [];
    runTimer = [];
    if nargin >= 3
        fileIdx = varargin{3};
    end
    if nargin >= 4
        totalFiles = varargin{4};
    end
    if nargin >= 5
        runTimer = varargin{5};
    end

    % 记录函数开始时间
    file_start_time = tic;
    
    % 【内存监控】显示初始内存状态
    try
        if ispc
            m = memory;
            fprintf('初始可用内存: %.1f GB\n', m.MemAvailableAllArrays/1024^3);
        end
    catch
        % 忽略内存监控错误
    end

    %% (2) 初始化参数
    % 加载参数
    params = config_params();
    fprintf('加载参数后: params.processing.max_frames = %d\n', params.processing.max_frames);
    % 是否由本函数启动并行池的标志（默认为 false）
    poolStartedHere = false;

    % 如果输出路径名遵循 'dopu_<N>layer[_<sigma>_<size>]' 约定，则从路径解析并覆盖部分参数
    try
            % 允许传入 run_output_base（例如 .../run_YYYY），若发现目录名以 'run_' 开头，则改用父目录名进行解析
            root_output_for_parsing = output_base;
            [~, possible_name, ~] = fileparts(root_output_for_parsing);
            if startsWith(possible_name, 'run_')
                root_output_for_parsing = fileparts(root_output_for_parsing);
            end
            [~, outbase_name, ~] = fileparts(root_output_for_parsing);
        
        % 根据输出目录名决定滤波模式：包含 'dopu' 则使用自适应DOPU (wovWinF = 0)，否则（例如 DDG）切换到固定高斯 (wovWinF = 1)
        if contains(outbase_name, 'dopu', 'IgnoreCase', true)
            params.mode.wovWinF = 0; % 自适应DOPU
            updates.mode.wovWinF = 0;
            fprintf('输出目录包含 "dopu"，设置 params.mode.wovWinF = 0 (自适应DOPU)\n');
            wrote_any = true;
        else
            params.mode.wovWinF = 1; % 固定高斯 / DDG 场景
            updates.mode.wovWinF = 1;
            fprintf('输出目录不包含 "dopu"，设置 params.mode.wovWinF = 1 (固定高斯 / DDG)\n');
            wrote_any = true;
        end
        
        % 先尝试匹配以 'dopu_' 开头的格式，再退回到通用 '<N>layer[_sigma_size]' 格式
        tok = regexp(outbase_name, 'dopu_(\d+)layer(?:_([\d\.]+)_([\d]+))?', 'tokens', 'once');
        pattern_used = '';
        if ~isempty(tok)
            t = tok; pattern_used = 'dopu_(N)layer[_sigma_size]';
        else
            tok = regexp(outbase_name, '(\d+)layer(?:_([\d\.]+)_([\d]+))?', 'tokens', 'once');
            if ~isempty(tok)
                t = tok; pattern_used = '(N)layer[_sigma_size]';
            else
                t = {};
            end
        end
        if ~isempty(t)
            fprintf('解析到 tokens (%s): %s\n', pattern_used, strjoin(t, ', '));
            % Avnum
            if numel(t) >= 1 && ~isempty(strtrim(t{1}))
                newAv = str2double(strtrim(t{1}));
                params.polarization.Avnum = newAv;
                updates.polarization.Avnum = newAv;
                fprintf('解析输出目录名 -> Avnum = %d\n', params.polarization.Avnum);
                wrote_any = true;
            end
            % h2 sigma and size
            if numel(t) >= 3 && ~isempty(strtrim(t{2})) && ~isempty(strtrim(t{3}))
                newSigma = str2double(strtrim(t{2}));
                sz = str2double(strtrim(t{3}));
                params.filters.h2_sigma = newSigma;
                params.filters.h2_size = [sz sz];
                % 重新生成高斯核
                params.filters.h2 = fspecial('gaussian', params.filters.h2_size, params.filters.h2_sigma);
                updates.filters.h2_sigma = params.filters.h2_sigma;
                updates.filters.h2_size = params.filters.h2_size;
                fprintf('解析输出目录名 -> h2_sigma = %g, h2_size = [%d %d]\n', params.filters.h2_sigma, sz, sz);
                wrote_any = true;
            end

            % 后备解析: 如果没有解析到 h2_sigma/h2_size，则尝试按下划线分割查找纯数字字段
            if (~isfield(updates,'filters') || ~isfield(updates.filters,'h2_sigma'))
                parts = strsplit(outbase_name, '_');
                numtokens = {};
                for pi = 1:numel(parts)
                    % 跳过包含 'layer' 的部分
                    if isempty(regexp(parts{pi}, '\d+layer', 'once'))
                        if ~isempty(regexp(parts{pi}, '^\d+(\.\d+)?$', 'once'))
                            numtokens{end+1} = parts{pi};
                        end
                    end
                end
                if numel(numtokens) >= 2
                    % 第一个数字作为 sigma，第二个数字作为 size
                    newSigma = str2double(numtokens{1});
                    sz = str2double(numtokens{2});
                    params.filters.h2_sigma = newSigma;
                    params.filters.h2_size = [sz sz];
                    params.filters.h2 = fspecial('gaussian', params.filters.h2_size, params.filters.h2_sigma);
                    updates.filters.h2_sigma = params.filters.h2_sigma;
                    updates.filters.h2_size = params.filters.h2_size;
                    fprintf('后备解析输出目录 -> h2_sigma = %g, h2_size = [%d %d]\n', newSigma, sz, sz);
                    wrote_any = true;
                end
            end

            % 如果解析到任意要写回的配置，则尝试写回并验证
            if wrote_any
                try
                    cfg_file = fullfile(fileparts(mfilename('fullpath')), 'function', 'config_params.m');
                    overwrite_config_params(updates, cfg_file);

                    % 读取写回后的文件并校验（使用简单子字符串匹配以避免正则转义问题）
                    cfg_text = fileread(cfg_file);
                    
                    % 验证 wovWinF
                    if isfield(updates, 'mode') && isfield(updates.mode, 'wovWinF')
                        expected_wovWinF = sprintf('params.mode.wovWinF = %d;', updates.mode.wovWinF);
                        if contains(cfg_text, expected_wovWinF)
                            fprintf('已确认写回: %s\n', expected_wovWinF);
                        else
                            warning('写回验证失败: %s 未在 %s 中更新', expected_wovWinF, cfg_file);
                        end
                    end
                    
                    % 验证 Avnum
                    if isfield(updates, 'polarization') && isfield(updates.polarization, 'Avnum')
                        expected = sprintf('params.polarization.Avnum = %d;', updates.polarization.Avnum);
                        if contains(cfg_text, expected)
                            fprintf('已确认写回: %s\n', expected);
                        else
                            warning('写回验证失败: %s 未在 %s 中更新', expected, cfg_file);
                        end
                    end
                    % 验证 h2_sigma 和 h2_size
                    if isfield(updates, 'filters') && isfield(updates.filters, 'h2_sigma')
                        expected_s = sprintf('params.filters.h2_sigma = %g;', updates.filters.h2_sigma);
                        if contains(cfg_text, expected_s)
                            fprintf('已确认写回: %s\n', expected_s);
                        else
                            warning('写回验证失败: %s 未在 %s 中更新', expected_s, cfg_file);
                        end
                    end
                    if isfield(updates, 'filters') && isfield(updates.filters, 'h2_size')
                        expected_sz = sprintf('params.filters.h2_size = [%d %d];', updates.filters.h2_size(1), updates.filters.h2_size(2));
                        if contains(cfg_text, expected_sz)
                            fprintf('已确认写回: %s\n', expected_sz);
                        else
                            warning('写回验证失败: %s 未在 %s 中更新', expected_sz, cfg_file);
                        end
                    end
                catch ME
                    warning('写回配置至 config_params.m 失败: %s', ME.message);
                end
            else
                fprintf('没有从输出目录解析到需要写回的配置项，跳过写回。\n');
            end
        else
            % 即使没有解析到 layer 参数，也要写回 wovWinF 配置
            if wrote_any
                try
                    cfg_file = fullfile(fileparts(mfilename('fullpath')), 'function', 'config_params.m');
                    overwrite_config_params(updates, cfg_file);
                    
                    % 验证 wovWinF 写回
                    cfg_text = fileread(cfg_file);
                    if isfield(updates, 'mode') && isfield(updates.mode, 'wovWinF')
                        expected_wovWinF = sprintf('params.mode.wovWinF = %d;', updates.mode.wovWinF);
                        if contains(cfg_text, expected_wovWinF)
                            fprintf('已确认写回: %s\n', expected_wovWinF);
                        else
                            warning('写回验证失败: %s 未在 %s 中更新', expected_wovWinF, cfg_file);
                        end
                    end
                catch ME
                    warning('写回配置至 config_params.m 失败: %s', ME.message);
                end
            end
        end
    catch ME
        warning('解析 output_base 名称失败: %s', ME.message);
    end

    % 启动并行池（安全版）
    % - 不再硬编码 worker 数
    % - 根据 local profile 的 NumWorkers 与系统物理核数选择
    % - 支持可选参数 params.parallel.maxWorkers（若在 config 中设置）
    if isempty(gcp('nocreate'))
        try
            c_local = parcluster('local');
            sysCores = feature('numcores');
            maxAllowed = min(c_local.NumWorkers, sysCores);
            % 如果在 config 中提供了并行最大 worker 数则使用它
            if isfield(params, 'parallel') && isfield(params.parallel, 'maxWorkers') && params.parallel.maxWorkers > 0
                targetWorkers = min(maxAllowed, params.parallel.maxWorkers);
                fprintf('配置中设置的最大worker数 (params.parallel.maxWorkers): %d\n', params.parallel.maxWorkers);
            else
                % 默认不要超过 8 个 worker，也不要超过可用核数
                targetWorkers = max(1, min(maxAllowed, 8));
                fprintf('使用默认worker数 (min(maxAllowed, 8)): %d\n', targetWorkers);
            end
            fprintf('尝试启动并行池 (local)，使用 %d 个 worker...\n', targetWorkers);
            parpool(c_local, targetWorkers);
            % 标记为本函数已启动并行池
            poolStartedHere = true;
        catch ME
            if isfield(ME, 'identifier') && ~isempty(ME.identifier)
                warning(ME.identifier, '%s\n将改为串行执行。', ME.message);
            else
                warning('parpool:startupFailed', '%s\n将改为串行执行。', ME.message);
            end
        end
    end


    % 直接用params.xxx结构体成员，无需单独赋值
    distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing

    %% （4）处理单个文件
    filename = input_file_path;
    [filepath, name, ext] = fileparts(filename);
    display_name = [name ext];

    fprintf('\n==================================================\n');
    fprintf('正在处理文件: %s\n', display_name);
    fprintf('文件路径: %s\n', filename);
    fprintf('==================================================\n');

    % This reads the parameter used for the data acquisition from *.oct* file
    fid               = fopen(filename);
    bob               = fread(fid, 1, 'uint32');
    SPL               = fread(fid, 1, 'double');
    nX                = fread(fid, 1, 'uint32'); %% number of Alines
    nY                = fread(fid, 1, 'uint32'); %% number of B-scans
    Boffset           = fread(fid, 1, 'uint32');
    Blength           = fread(fid, 1, 'uint32') + 1;
    Xcenter           = fread(fid, 1, 'double');
    Xspan             = fread(fid, 1, 'double');
    Ycenter           = fread(fid, 1, 'double');
    Yspan             = fread(fid, 1, 'double');
    frame_per_pos     = fread(fid, 1, 'uint32'); %% repetition of B-scans
    n_dataset         = fread(fid, 1, 'uint32'); %% repetion of volume scan
    ProtMode          = fread(fid, 1, 'uint32');
    fseek(fid, 4, 'cof');%v10
    sizeBck           = fread(fid, 1, 'uint32');
    Bck1              = fread(fid, sizeBck, 'int16');
    sizeKES           = fread(fid, 2, 'uint32');
    KES1              = (fread(fid, sizeKES(2), 'double'))' * sizeKES(2);
    sizeBck           = fread(fid, 1, 'uint32');
    Bck2              = fread(fid, sizeBck, 'int16');
    sizeKES           = fread(fid, 2, 'uint32'); %new 4096 cal
    KES2              = (fread(fid, sizeKES(2), 'double'))' * sizeKES(2);
    disp_coef         = fread(fid, 1, 'double'); %%dispersion coefficient
    nR                = frame_per_pos;
    IMGheight         = floor(Blength / 2);
    kmat              = linspace(-0.5, 0.5, Blength)'.^2;
    phV               = exp(1i .* (kmat .* disp_coef));
    nY                = floor(nY / nR);
    K                 = 1:Blength;
    use_autoRg        = 1;
    RgFlow            = [60 110];
    jsurf             = zeros(4);
    IMcropRg          = 1:(SPL/2);
    nZcrop            = numel(IMcropRg);
    imshowrgZ         = 1:nZcrop;
    jusamp            = zeros(nZcrop, nX, 4);
    winG              = tukeywin(Blength, 0.25);
    %(a)subWins: params for split spectrum DOPU
    nWin              = 9;
    winL              = 2 * Blength / (nWin + 1);
    winG              = tukeywin(winL, 0.25);
    winG_whole        = tukeywin(Blength, 0.25); % window for whole spectrum
    windex            = 1 : winL / 2 : Blength;

    % 添加帧数限制功能
    if params.processing.max_frames > 0
        original_nY = nY;
        nY = min(nY, params.processing.max_frames);
        fprintf('帧数限制: %d -> %d (最大允许: %d)\n', original_nY, nY, params.processing.max_frames);
    end

    if params.processing.useref == 1 %use the first 50k a-lines to calc ref
        fseek(fid, bob, 'bof');
        n_ref = floor(min(50000, floor(nY * nX / 2) * 2) / nX);
        Ref_ch1 = zeros(SPL, nX, n_ref);
        Ref_ch2 = zeros(SPL, nX, n_ref);
        for i_ref = 1:n_ref
            fseek(fid, 4, 'cof');
            BB = reshape(fread(fid, SPL * nX * 2, 'int16'), [SPL, nX * 2]);
            Ref_ch1(:, :, i_ref) = BB(:, 1:nX);
            Ref_ch2(:, :, i_ref) = BB(:, nX + 1:end);
        end
        Ref_ch1 = repmat(mean(mean(Ref_ch1, 2), 3), [1, nX]);
        Ref_ch2 = repmat(mean(mean(Ref_ch2, 2), 3), [1, nX]);
    elseif params.processing.useref == -1
        Ref_ch1 = 0;
        Ref_ch2 = 0;
    else
        Ref_ch1 = repmat(Bck1, [1, nX]);
        Ref_ch2 = repmat(Bck2, [1, nX]);
    end


    %(b) Calculates the Cumulative Stokes parameters I,Q,U,V
    for nr = nR
        % 生成清晰的输出目录（每个输入文件单独子目录）
        % 默认行为： output_base/<filename>；若存在则在名称后附加时间戳，除非在配置中设置 params.processing.overwrite_output = true
        try
            if ~isfield(params, 'processing') || ~isfield(params.processing, 'output_subdir_unique')
                % 默认启用每文件子目录
                params.processing.output_subdir_unique = true;
            end
            if ~isfield(params.processing, 'overwrite_output')
                params.processing.overwrite_output = false;
            end
        catch
            % 忽略错误，使用默认行为
        end

        if params.processing.output_subdir_unique
            base_name = name;
            candidate_dir = fullfile(output_base, base_name);
            if exist(candidate_dir, 'dir')
                if params.processing.overwrite_output
                    foutputdir = candidate_dir;
                    fprintf('输出目录已存在，且配置允许覆盖：%s\n', foutputdir);
                else
                    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                    foutputdir = fullfile(output_base, [base_name '_' timestamp]);
                    fprintf('输出目录已存在，使用新的目录: %s\n', foutputdir);
                end
            else
                foutputdir = candidate_dir;
            end
            if ~exist(foutputdir, 'dir')
                mkdir(foutputdir);
            end
        else
            % 使用共享 output_base 目录
            foutputdir = output_base;
            if ~exist(foutputdir, 'dir'), mkdir(foutputdir); end
        end

        czrg = 1:320; % set z range
        topLines = ones(nX, nY);
        if params.processing.hasSeg
            matFilePath = fullfile(filepath, [name, '.mat']);
            if exist(matFilePath, 'file')
                try
                    load(matFilePath, 'topLines');
                    topLines = double(topLines);
                    
                    % ========== 自动修复 topLines 维度 ==========
                    [loaded_dim1, loaded_dim2] = size(topLines);
                    
                    % 检测是否需要转置
                    if loaded_dim1 ~= nX && loaded_dim2 == nX
                        topLines = topLines';
                        [loaded_dim1, loaded_dim2] = size(topLines);
                    end
                    
                    % 按实际帧数切片（处理 max_frames 限制）
                    if loaded_dim1 == nX && loaded_dim2 > nY
                        topLines = topLines(:, 1:nY);
                    elseif loaded_dim1 == nX && loaded_dim2 < nY
                        topLines_padded = ones(nX, nY);
                        topLines_padded(:, 1:loaded_dim2) = topLines;
                        topLines = topLines_padded;
                    elseif loaded_dim1 ~= nX
                        warning('topLines 维度不匹配，将重置为默认值');
                        topLines = ones(nX, nY);
                        params.processing.hasSeg = 0;
                    end
                    % ==========================================================
                    
                    % 注意: 无效边界值(0/负数/NaN)将在后续统一修复为 nZ (底部)
                    fprintf('成功加载分割结果文件: %s\n', matFilePath);
                catch ME
                    warning('加载分割结果文件失败: %s\n错误: %s\n将自动进行表面分割', matFilePath, ME.message);
                    params.processing.hasSeg = 0;  % 自动切换到不使用预分割模式
                end
            else
                warning('分割结果文件不存在: %s\n将自动进行表面分割', matFilePath);
                params.processing.hasSeg = 0;  % 自动切换到不使用预分割模式
            end
        end
        if params.range.setZrg
            czrg = czrg(1:round(params.range.setZrg));
        end
        nZcrop = numel(czrg);
        nZ = nZcrop;
        
        % ========== 修复无效 Topline 边界值 ==========
        if exist('topLines', 'var')
            invalid_mask = (topLines <= 0) | isnan(topLines);
            if any(invalid_mask(:))
                topLines(invalid_mask) = 1;
                fprintf('已修复 %d 个无效 topLines 值 (<=0 或 NaN)\n', sum(invalid_mask(:)));
            end
        end
        % =========================================================
        
        % (c) creat array to store results
        LA_c_cfg1_avg = zeros(nZcrop - 20, nX, 3, nY);
        PhR_c_cfg1_avg = zeros(nZcrop - 20, nX, nY);
        cumLA_cfg1_avg = zeros(nZcrop - 20, nX, 3, nY);
        LA_c_cfg1_eig = zeros(nZcrop - 20, nX, 3, nY);
        % PhR_c_cfg1_eig = zeros(nZcrop - 20, nX, nY);
        % cumLA_cfg1_eig = zeros(nZcrop - 20, nX, 3, nY);

        LA_Ms_cfg1_rmBG = LA_c_cfg1_avg;
        PhR_Ms_cfg1_rmBG = cumLA_cfg1_avg;
        cumLA_Ms_cfg1_rmBG = LA_c_cfg1_eig;
        Smap_avg = zeros(numel(czrg), nX, 3, nY);
        Smap_rep1 = zeros(numel(czrg), nX, 3, nY);
        Strus = zeros(numel(czrg), nX, nY);
        Stru_OAC = Strus;
        dopu_splitSpectrum = zeros(numel(czrg), nX, nY);
        
        % ========== 加载后巩膜边界数据（支持多文件自动匹配）==========
        % 加载逻辑优先级：
        %   1. 首先尝试从输入文件夹下的 scaler_mat 子文件夹中查找匹配的 mat 文件
        %   2. 如果 scaler_mat 中找不到，则回退到配置文件中指定的路径
        sclera_boundary_data = [];  % 初始化为空
        has_sclera_boundary = false;
        sclera_boundary_file = '';
        
        % 检查是否启用后巩膜边界处理
        use_sclera = true;  % 默认启用
        if isfield(params, 'files') && isfield(params.files, 'use_sclera_boundary')
            use_sclera = params.files.use_sclera_boundary;
        end
        
        if ~use_sclera
            fprintf('\n后巩膜边界处理已禁用 (params.files.use_sclera_boundary = 0)\n');
        else
            % 获取 scaler_mat 文件夹名称（可配置）
            scaler_folder_name = 'scaler_mat';  % 默认名称
            if isfield(params, 'files') && isfield(params.files, 'sclera_mat_folder') && ~isempty(params.files.sclera_mat_folder)
                scaler_folder_name = params.files.sclera_mat_folder;
            end
            
            % 方式1: 自动从 scaler_mat 文件夹中查找匹配的后巩膜边界文件
            scaler_mat_folder = fullfile(filepath, scaler_folder_name);
            if exist(scaler_mat_folder, 'dir')
                fprintf('\n========== 搜索后巩膜边界数据 (%s) ==========\n', scaler_folder_name);
                fprintf('搜索目录: %s\n', scaler_mat_folder);
                
                % 获取当前处理文件的基础名（不含扩展名），用于匹配
                [~, current_base_name, ~] = fileparts(name);
                
                % 搜索所有 .mat 文件
                mat_files_in_scaler = dir(fullfile(scaler_mat_folder, '*.mat'));
                
                if ~isempty(mat_files_in_scaler)
                    fprintf('找到 %d 个 MAT 文件\n', length(mat_files_in_scaler));
                    
                    % 尝试多种匹配策略
                    matched_file = '';
                    
                    % 策略1: 精确匹配（文件名相同）
                    for k = 1:length(mat_files_in_scaler)
                        [~, mat_base_name, ~] = fileparts(mat_files_in_scaler(k).name);
                        if strcmpi(mat_base_name, current_base_name)
                            matched_file = fullfile(scaler_mat_folder, mat_files_in_scaler(k).name);
                            fprintf('精确匹配成功: %s\n', mat_files_in_scaler(k).name);
                            break;
                        end
                    end
                    
                    % 策略2: 模糊匹配（文件名包含当前处理文件的基础名）
                    if isempty(matched_file)
                        for k = 1:length(mat_files_in_scaler)
                            [~, mat_base_name, ~] = fileparts(mat_files_in_scaler(k).name);
                            if contains(mat_base_name, current_base_name, 'IgnoreCase', true) || ...
                               contains(current_base_name, mat_base_name, 'IgnoreCase', true)
                                matched_file = fullfile(scaler_mat_folder, mat_files_in_scaler(k).name);
                                fprintf('模糊匹配成功: %s\n', mat_files_in_scaler(k).name);
                                break;
                            end
                        end
                    end
                    
                    % 策略3: 如果只有一个文件，直接使用
                    if isempty(matched_file) && length(mat_files_in_scaler) == 1
                        matched_file = fullfile(scaler_mat_folder, mat_files_in_scaler(1).name);
                        fprintf('仅有一个文件，自动使用: %s\n', mat_files_in_scaler(1).name);
                    end
                    
                    % 策略4: 按文件索引匹配（如果有多个文件且有 fileIdx）
                    if isempty(matched_file) && ~isempty(fileIdx) && fileIdx <= length(mat_files_in_scaler)
                        matched_file = fullfile(scaler_mat_folder, mat_files_in_scaler(fileIdx).name);
                        fprintf('按索引匹配 (第 %d 个文件): %s\n', fileIdx, mat_files_in_scaler(fileIdx).name);
                    end
                    
                    if ~isempty(matched_file)
                        sclera_boundary_file = matched_file;
                    else
                        fprintf('未能在 %s 中找到匹配的后巩膜边界文件\n', scaler_folder_name);
                        % 列出可用文件供参考
                        fprintf('可用文件列表:\n');
                        for k = 1:min(5, length(mat_files_in_scaler))
                            fprintf('  [%d] %s\n', k, mat_files_in_scaler(k).name);
                        end
                        if length(mat_files_in_scaler) > 5
                            fprintf('  ... 等共 %d 个文件\n', length(mat_files_in_scaler));
                        end
                    end
                else
                    fprintf('%s 文件夹中没有找到 MAT 文件\n', scaler_folder_name);
                end
            end
            
            % 方式2: 如果 scaler_mat 中未找到，回退到配置文件中指定的路径
            if isempty(sclera_boundary_file) && isfield(params, 'files') && isfield(params.files, 'sclera_boundary_path') && ~isempty(params.files.sclera_boundary_path)
                fprintf('回退到配置文件中指定的后巩膜边界路径\n');
                sclera_boundary_file = params.files.sclera_boundary_path;
            end
        end
        
        % 加载后巩膜边界数据
        if use_sclera && ~isempty(sclera_boundary_file) && exist(sclera_boundary_file, 'file')
            fprintf('\n========== 加载后巩膜边界数据 ==========\n');
            fprintf('文件路径: %s\n', sclera_boundary_file);
            
            try
                % 确定变量名（默认与topLines相同的变量名，或使用配置的变量名）
                if isfield(params, 'files') && isfield(params.files, 'sclera_boundary_var') && ~isempty(params.files.sclera_boundary_var)
                    var_name = params.files.sclera_boundary_var;
                else
                    % 默认尝试常见的变量名
                    var_name = 'bottomLines';  % 默认名称
                    % 检查文件中是否存在该变量，如果不存在则尝试其他名称
                    mat_info = whos('-file', sclera_boundary_file);
                    var_names = {mat_info.name};
                    
                    common_names = {'bottomLines', 'scleraLines', 'boundary', 'layer', 'sclera_line', 'topLines'};
                    found = false;
                    for k = 1:length(common_names)
                        if ismember(common_names{k}, var_names)
                            var_name = common_names{k};
                            found = true;
                            break;
                        end
                    end
                    
                    % 如果都没找到，使用第一个变量
                    if ~found && ~isempty(var_names)
                        var_name = var_names{1};
                    end
                end
                
                % 加载边界数据（参考topLines的加载方式）
                loaded_data = load(sclera_boundary_file, var_name);
                
                % 直接获取加载的变量
                if isfield(loaded_data, var_name)
                    sclera_boundary_data = loaded_data.(var_name);
                    fprintf('成功从变量 "%s" 加载边界数据\n', var_name);
                    
                    % 转换为double（与topLines处理一致）
                    sclera_boundary_data = double(sclera_boundary_data);
                    
                    % 注意: 无效边界值将在后续统一修复
                    
                    % 显示加载信息
                    fprintf('边界数据维度: [%s]\n', num2str(size(sclera_boundary_data)));
                    fprintf('边界范围: %.1f ~ %.1f\n', min(sclera_boundary_data(:)), max(sclera_boundary_data(:)));
                    fprintf('后巩膜边界数据加载成功！\n');
                    fprintf('======================================\n\n');
                    
                    has_sclera_boundary = true;
                else
                    warning('变量 "%s" 加载失败', var_name);
                    sclera_boundary_data = [];
                    has_sclera_boundary = false;
                end
                
            catch ME
                warning('加载后巩膜边界数据失败: %s', ME.message);
                fprintf('错误详情: %s\n', ME.getReport());
                sclera_boundary_data = [];
                has_sclera_boundary = false;
            end
        elseif ~isempty(sclera_boundary_file)
            warning('后巩膜边界文件不存在: %s', sclera_boundary_file);
        else
            fprintf('未配置后巩膜边界数据路径，跳过加载\n');
        end
        
        % ========== 修复无效后巩膜边界值 ==========
        if has_sclera_boundary && ~isempty(sclera_boundary_data)
            invalid_mask_sclera = (sclera_boundary_data <= 0) | isnan(sclera_boundary_data);
            if any(invalid_mask_sclera(:))
                sclera_boundary_data(invalid_mask_sclera) = 1;
                fprintf('已修复 %d 个无效后巩膜边界值\n', sum(invalid_mask_sclera(:)));
            end
        end
        % ====================================================================
        
        fprintf('开始处理 %d 个B-Scan...\n', nY);
        % ========================= 分批流式处理架构 =========================
        % 修改目标:
        %   1. 移除全量预加载（all_Bs1/all_Bs2），改为分批处理
        %   2. 数据类型从 double 改为 single（内存减半）
        %   3. 串行读取 + 并行计算模式，避免 I/O 竞争
        % ====================================================================
        
        % 配置日志等级（0=静默,1=简短,2=详尽）
        if exist('verbose', 'var') == 0
            if isfield(params, 'logging') && isfield(params.logging, 'verbose')
                verbose = params.logging.verbose;
            elseif isfield(params, 'processing') && isfield(params.processing, 'verbose')
                verbose = params.processing.verbose;
            else
                verbose = 0; % 默认静默
            end
        end
        
        % 设置批处理大小（可配置）
        if isfield(params, 'parallel') && isfield(params.parallel, 'batchSize')
            BatchSize = params.parallel.batchSize;
        else
            BatchSize = 5; % 默认每批处理 50 帧
        end
        
        fprintf('\n========== 分批流式处理模式 ==========\n');
        fprintf('总帧数: %d, 批大小: %d, 批次数: %d\n', nY, BatchSize, ceil(nY/BatchSize));
        fprintf('使用 single 精度（内存优化）\n');
        fprintf('=====================================\n\n');
        
        % 记录总体开始时间
        total_batch_start_time = tic;
        total_batches = ceil(nY/BatchSize);
        
        % 分批处理主循环
        for batch_start = 1:BatchSize:nY
            batch_end = min(batch_start + BatchSize - 1, nY);
            batch_size_actual = batch_end - batch_start + 1;
            current_batch = ceil(batch_start/BatchSize);
            
            fprintf('\n--- 批次 %d/%d: 处理帧 %d-%d (%d 帧) ---\n', ...
                current_batch, total_batches, ...
                batch_start, batch_end, batch_size_actual);
            
            % ========== 步骤1: 串行读取当前批次数据（使用 single 精度）==========
            t_load = tic;
            fid_read = fopen(filename);
            
            % 预分配当前批次的内存（single 类型，内存减半）
            batch_Bs1 = zeros(Blength, nX, nr, batch_size_actual, 'single');
            batch_Bs2 = zeros(Blength, nX, nr, batch_size_actual, 'single');
            
            fprintf('  读取数据: ');
            for iY_local = 1:batch_size_actual
                % 显示读取进度（每10%更新一次）
                if mod(iY_local, max(1, floor(batch_size_actual/5))) == 0 || iY_local == batch_size_actual
                    print_progress(iY_local, batch_size_actual, '读取', 20);
                end
                
                iY_global = batch_start + iY_local - 1;
                fseek(fid_read, bob + (SPL * nX * 2 + 2) * 2 * (nR * (iY_global - 1)), 'bof');
                for Ic = 1:nr
                    fseek(fid_read, 4, 'cof');
                    B = single(reshape(fread(fid_read, SPL * nX * 2, 'int16'), [SPL, nX * 2]));
                    batch_Bs1(:, :, Ic, iY_local) = (B(:, 1:nX) - single(Ref_ch1));
                    batch_Bs2(:, :, Ic, iY_local) = (B(:, nX + 1:end) - single(Ref_ch2));
                end
            end
            fclose(fid_read);
            loadTime = toc(t_load);
            fprintf('  数据读取完成，耗时: %.2f 秒\n', loadTime);
            
            % ========== 步骤2: 并行处理当前批次 ==========
            t_proc = tic;
            
            % 为当前批次预分配临时数组（parfor 要求使用循环变量作为索引）
            batch_dopu = zeros(numel(czrg), nX, batch_size_actual);
            batch_Strus = zeros(numel(czrg), nX, batch_size_actual);
            batch_Smap_avg = zeros(numel(czrg), nX, 3, batch_size_actual);
            batch_Smap_rep1 = zeros(numel(czrg), nX, 3, batch_size_actual);
            batch_topLines = topLines(:, batch_start:batch_end);  % 当前批次的 topLines
            batch_LA_c_cfg1_avg = zeros(nZcrop - 20, nX, 3, batch_size_actual);
            batch_PhR_c_cfg1_avg = zeros(nZcrop - 20, nX, batch_size_actual);
            batch_cumLA_cfg1_avg = zeros(nZcrop - 20, nX, 3, batch_size_actual);
            batch_LA_Ms_cfg1_rmBG = zeros(nZcrop - 20, nX, 3, batch_size_actual);
            batch_PhR_Ms_cfg1_rmBG = zeros(nZcrop - 20, nX, batch_size_actual);
            batch_cumLA_Ms_cfg1_rmBG = zeros(nZcrop - 20, nX, 3, batch_size_actual);
            
            % 将需要的参数提取为局部变量（避免 parfor 中的广播变量问题）
            local_hasSeg = params.processing.hasSeg;
            local_do_PhComp = params.processing.do_PhComp;
            local_do_combined = params.dopu.do_combined;
            local_do_spatial = params.dopu.do_spatial;
            local_kRL_cfg1 = params.polarization.kRL_cfg1;
            local_kRU_cfg1 = params.polarization.kRU_cfg1;
            local_h1 = params.filters.h1;
            local_h2 = params.filters.h2;
            local_Avnum = params.polarization.Avnum;
            local_wovWinF = params.mode.wovWinF;
            local_enableDopuPhaseSupp = params.polarization.enableDopuPhaseSupp;
            local_enableOutputAdaptive = params.filters.enable_output_adaptive;
            local_kRL_output = params.filters.kRL_output;
            local_kRU_output = params.filters.kRU_output;
            local_outputDopuThreshold = params.filters.output_dopu_threshold;
            local_enableBottomLayerPhaseReduction = params.filters.enable_bottom_layer_phase_reduction;
            local_bottomLayerDepth = params.filters.bottom_layer_depth;
            local_adaptiveFilterBottomDepth = params.filters.adaptive_filter_bottom_depth;
            % 从配置中读取底层DOPU阈值（用于Step 3.6），如果不存在则使用默认0.5
            if isfield(params.filters, 'bottom_dopu_threshold')
                local_bottomDopuThreshold = params.filters.bottom_dopu_threshold;
            else
                local_bottomDopuThreshold = 0.5;
            end
            local_phV = single(phV);
            local_winG = single(winG);
            local_winG_whole = single(winG_whole);
            
            parfor iY_local = 1:batch_size_actual
                Bs1 = squeeze(batch_Bs1(:, :, :, iY_local));
                Bs2 = squeeze(batch_Bs2(:, :, :, iY_local));
                if local_do_PhComp == 1
                    Bd1 = real(hilbert(Bs1) .* local_phV);  % dispersion correction (转换为single)
                    Bd2 = real(hilbert(Bs2) .* local_phV);
                else
                    Bd1 = Bs1;
                    Bd2 = Bs2;
                end

                %% split spectrum / spatial / combined DOPU
                if local_do_combined
                    % 使用组合DOPU (分裂谱 + 空间)
                    dopu_result = calculateCombinedDOPU(...
                        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, local_winG, czrg);
                    batch_dopu(:, :, iY_local) = dopu_result;
                    dopu_ss = dopu_result;
                elseif local_do_spatial
                    % 使用空间DOPU
                    Bimg1_wholeStr = fft(Bd1 .* local_winG_whole, SPL, 1);
                    Bimg2_wholeStr = fft(Bd2 .* local_winG_whole, SPL, 1);
                    IMG1_wholeStr = Bimg1_wholeStr(czrg, :, :);
                    IMG2_wholeStr = Bimg2_wholeStr(czrg, :, :);
                    dopu_result = calculateSpatialDOPU(IMG1_wholeStr, IMG2_wholeStr, params);
                    batch_dopu(:, :, iY_local) = dopu_result;
                    dopu_ss = dopu_result;
                else
                    % 使用分裂谱DOPU (默认)
                    [dopu_result, dopu_ss] = calculateSplitSpectrumDOPU(...
                        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, local_winG, czrg);
                    batch_dopu(:, :, iY_local) = dopu_result;
                end

                %% whole spectrum nZ*nX*nR==> fft(complex)
                Bimg1_wholeStr = fft(Bd1 .* local_winG_whole, SPL, 1);
                Bimg2_wholeStr = fft(Bd2 .* local_winG_whole, SPL, 1);
                IMG1_wholeStr = Bimg1_wholeStr(czrg, :, :);
                IMG2_wholeStr = Bimg2_wholeStr(czrg, :, :);
 
                %% Struc Stokes, and OAC
                [wS0, wS1, wS2, wS3] = cumulativeQUV(IMG1_wholeStr, IMG2_wholeStr);
                wQ = wS1 ./ wS0;
                wU = wS2 ./ wS0;
                wV = wS3 ./ wS0;
                strLin = mean(wS0, 3);
                batch_Strus(:, :, iY_local) = 20 * log10(strLin);
                batch_Smap_avg(:, :, :, iY_local) = cat(3, mean(wQ, 3), mean(wU, 3), mean(wV, 3));
                batch_Smap_rep1(:, :, :, iY_local) = cat(3, wQ(:, :, 1), wU(:, :, 1), wV(:, :, 1));
                
                % 表面分割（如果没有预分割）
                local_topLine = batch_topLines(:, iY_local);
                if ~local_hasSeg
                    strOAC = calOAC(strLin); 
                    local_topLine = surf_seg(strOAC, 0.25) + 2;
                    batch_topLines(:, iY_local) = local_topLine;
                end
                
                %% drLA and drPhR
                dopu_splitSpec_M = squeeze(mean(dopu_ss, 3));
                IMG_ch1 = squeeze(mean(IMG1_wholeStr, 3));
                IMG_ch2 = squeeze(mean(IMG2_wholeStr, 3));
                [LA_tmp, PhR_tmp, cumLA_tmp, LA_rmBG_tmp, PhR_rmBG_tmp, cumLA_rmBG_tmp] = ...
                    calLAPhRALL(IMG_ch1, IMG_ch2, local_topLine, dopu_splitSpec_M, ...
                        local_kRL_cfg1, local_kRU_cfg1, local_h1, local_h2, ...
                        local_Avnum, local_wovWinF, local_enableDopuPhaseSupp, verbose, ...
                        local_enableOutputAdaptive, local_kRL_output, local_kRU_output, ...
                        local_outputDopuThreshold, local_enableBottomLayerPhaseReduction, local_bottomLayerDepth, local_adaptiveFilterBottomDepth, local_bottomDopuThreshold);
                
                % 存储到批次临时数组（使用 iY_local 作为索引）
                batch_LA_c_cfg1_avg(:, :, :, iY_local) = LA_tmp;
                batch_PhR_c_cfg1_avg(:, :, iY_local) = PhR_tmp;
                batch_cumLA_cfg1_avg(:, :, :, iY_local) = cumLA_tmp;
                batch_LA_Ms_cfg1_rmBG(:, :, :, iY_local) = LA_rmBG_tmp;
                batch_PhR_Ms_cfg1_rmBG(:, :, iY_local) = PhR_rmBG_tmp;
                batch_cumLA_Ms_cfg1_rmBG(:, :, :, iY_local) = cumLA_rmBG_tmp;
            end % parfor end（当前批次处理完成）
            
            % ========== 步骤3: 将批次结果复制到全局数组 ==========
            dopu_splitSpectrum(:, :, batch_start:batch_end) = batch_dopu;
            Strus(:, :, batch_start:batch_end) = batch_Strus;
            Smap_avg(:, :, :, batch_start:batch_end) = batch_Smap_avg;
            Smap_rep1(:, :, :, batch_start:batch_end) = batch_Smap_rep1;
            topLines(:, batch_start:batch_end) = batch_topLines;
            LA_c_cfg1_avg(:, :, :, batch_start:batch_end) = batch_LA_c_cfg1_avg;
            PhR_c_cfg1_avg(:, :, batch_start:batch_end) = batch_PhR_c_cfg1_avg;
            cumLA_cfg1_avg(:, :, :, batch_start:batch_end) = batch_cumLA_cfg1_avg;
            LA_Ms_cfg1_rmBG(:, :, :, batch_start:batch_end) = batch_LA_Ms_cfg1_rmBG;
            PhR_Ms_cfg1_rmBG(:, :, batch_start:batch_end) = batch_PhR_Ms_cfg1_rmBG;
            cumLA_Ms_cfg1_rmBG(:, :, :, batch_start:batch_end) = batch_cumLA_Ms_cfg1_rmBG;
            
            procTime = toc(t_proc);
            fprintf('  并行计算完成，耗时: %.2f 秒\n', procTime);
            fprintf('  批次总耗时: %.2f 秒（读取 %.2f + 计算 %.2f）\n', loadTime + procTime, loadTime, procTime);
            
            % 显示总体进度
            total_elapsed = toc(total_batch_start_time);
            avg_time_per_batch = total_elapsed / current_batch;
            remaining_batches = total_batches - current_batch;
            eta_seconds = avg_time_per_batch * remaining_batches;
            fprintf('\n批次进度: ');
            print_progress(current_batch, total_batches, '', 30);
            fprintf('\n');
            
            % 同步更新全局（全部文件）进度（如果提供了 fileIdx/totalFiles/runTimer）
            if ~isempty(fileIdx) && ~isempty(totalFiles) && ~isempty(runTimer) && totalFiles > 1
                globalTotal = totalFiles * 100;
                % 在文件内按批次细分到 0-100 区间
                globalCurrent = (fileIdx - 1) * 100 + round(current_batch / total_batches * 100);
                elapsed_all = toc(runTimer);
                if globalCurrent > 0
                    avg_unit = elapsed_all / max(1, globalCurrent);
                    eta_all = avg_unit * max(0, globalTotal - globalCurrent);
                else
                    eta_all = NaN;
                end
                prefix_total = sprintf('总进度 %d/%d %s', fileIdx, totalFiles, name);
                fprintf('\n'); % 确保在新行显示总进度，不覆盖当前行输出
                print_progress(globalCurrent, globalTotal, prefix_total, 40, 'blocks', eta_all);
            end
            
            % ========== 步骤4: 清空批次数据（释放内存）==========
            clear batch_Bs1 batch_Bs2 batch_dopu batch_Strus batch_Smap_avg batch_Smap_rep1;
            clear batch_topLines batch_LA_c_cfg1_avg batch_PhR_c_cfg1_avg batch_cumLA_cfg1_avg;
            clear batch_LA_Ms_cfg1_rmBG batch_PhR_Ms_cfg1_rmBG batch_cumLA_Ms_cfg1_rmBG;
            
        end % 分批循环结束
        
        % 完成当前文件的所有批次后，确保更新全局进度为该文件完成
        if ~isempty(fileIdx) && ~isempty(totalFiles) && ~isempty(runTimer) && totalFiles > 1
            globalTotal = totalFiles * 100;
            globalCurrent = fileIdx * 100;
            elapsed_all = toc(runTimer);
            if globalCurrent > 0
                avg_unit = elapsed_all / max(1, globalCurrent);
                eta_all = avg_unit * max(0, globalTotal - globalCurrent);
            else
                eta_all = NaN;
            end
            prefix_total = sprintf('总进度 %d/%d %s', fileIdx, totalFiles, name);
            print_progress(globalCurrent, globalTotal, prefix_total, 40, 'blocks', eta_all);
        end
        fprintf('\n========== 所有批次处理完成！==========\n');
        fprintf('总耗时: %.2f 秒 (%.2f 分钟)\n', toc(total_batch_start_time), toc(total_batch_start_time)/60);
        fprintf('======================================\n\n');
    end
    fclose all;
    %% save results: strus(flow),stokes,oac
    if params.tiff.saveDicom
        % 创建 dcm 子文件夹
        dcm_dir = fullfile(foutputdir, 'dcm');
        if ~exist(dcm_dir, 'dir'), mkdir(dcm_dir); end
        
        % 创建低质量 dcm 子文件夹（用于快速浏览）- 根据配置决定是否启用
        if params.tiff.save_low_quality_dcm
            dcm_low_quality_dir = fullfile(foutputdir, 'dcm_low_quality');
            if ~exist(dcm_low_quality_dir, 'dir'), mkdir(dcm_low_quality_dir); end
            low_quality_scale = params.tiff.low_quality_scale;
        else
            dcm_low_quality_dir = '';  % 空字符串表示不生成低质量版本
            low_quality_scale = 0;
        end
        
        % 修复索引越界问题
        slice_index = min(100, floor(size(Strus,3)/2));  % 如果不足100层，则取中间层
        if slice_index < 1
            slice_index = 1;  % 至少取第一层
        end
        SS=Strus(:,:,slice_index);strUrg = max(SS(:))-5;strLrg = min(SS(:))+5;
        % 只处理前300层深度
        nZ_save = min(300, size(Strus, 1));
        for i = 1:size(Strus, 3)
            SS1 = Strus(1:nZ_save, :, i);  % 只取前300层
            Struc(1:nZ_save, :, :, i) = (SS1 - strLrg) ./ (strUrg - strLrg);
        end
        % 保存原始结构图像到 dcm 文件夹
        dicomwrite_with_lowquality(uint8(255 * (Struc)), fullfile(dcm_dir, [name, '_1-1_Struc.dcm']), dcm_low_quality_dir, low_quality_scale);

        % 创建带边界的结构图像副本（只保留前300层）
        Struc_with_boundary = Struc(1:nZ_save, :, :, :);
        % 对Struc应用边界处理
        for iY = 1:size(Struc_with_boundary, 4)
            for iX = 1:size(Struc_with_boundary, 2)
                surface_pos = round(topLines(iX, iY));
                
                if surface_pos > nZ_save
                    Struc_with_boundary(1:nZ_save, iX, :, iY) = 0;
                elseif surface_pos > 1
                    Struc_with_boundary(1:surface_pos-1, iX, :, iY) = 0;
                end
            end
        end
        dicomwrite_with_lowquality(uint8(255 * (Struc(1:nZ_save, :, :, :))), fullfile(dcm_dir, [name, '_1-1_Struc.dcm']), dcm_low_quality_dir, low_quality_scale);
        dicomwrite_with_lowquality(uint8(255 * (Struc_with_boundary)), fullfile(dcm_dir, [name, '_1-2_Struc_with_boundary.dcm']), dcm_low_quality_dir, low_quality_scale);
        % dicomwrite(uint8(255 * (Smap_rep1 / 2 + 0.5)), fullfile(dcm_dir, [name, '_1-3_1rep-Stokes.dcm']));
        dicomwrite_with_lowquality(uint8(255 * (Smap_avg(1:nZ_save, :, :, :) / 2 + 0.5)), fullfile(dcm_dir, [name, '_1-3_4rep-Stokes.dcm']), dcm_low_quality_dir, low_quality_scale);

        % 对1-4 DOPU图像应用阈值过滤（只保留前nZ_save层）
        dopu_thresholded = dopu_splitSpectrum(1:nZ_save, :, :);
        dopu_thresholded(dopu_thresholded <= 0.4) = 0;  % 小于等于0.4的设为0
        
        % 对DOPU图像应用表面边界mask（边界以上区域设为0）
        dopu_with_boundary = dopu_thresholded;
        for iY = 1:nY
            for iX = 1:size(dopu_with_boundary, 2)
                surface_pos = round(topLines(iX, iY));
                
                if surface_pos > nZ_save
                    dopu_with_boundary(1:nZ_save, iX, iY) = 0;
                elseif surface_pos > 1
                    dopu_with_boundary(1:surface_pos-1, iX, iY) = 0;
                end
            end
        end

        dicomwrite_with_lowquality(uint8(255 * (permute(dopu_with_boundary, [1 2 4 3]))), fullfile(dcm_dir, [name, '_1-4_dopu_SS.dcm']), dcm_low_quality_dir, low_quality_scale);
        
        if ~params.processing.hasSeg
            save(fullfile(filepath, [name, '.mat']), 'topLines', 'czrg');
        end
        rotAngle = 440;
        % 仅输出 cfg1 avg 相关结果（当前仅支持 avg 路径）
        PRRrg = params.polarization.PRRrg;
        writematrix(PRRrg, fullfile(dcm_dir, [name, '_2-0_PhRRg.txt']));
        
        % 生成彩色编码图像（直接写入 DCM，不生成多帧 TIFF）
        fprintf('生成彩色编码图像...\n');
        
        % 预分配数组
        cumLA_cfg_hsv = zeros(size(cumLA_cfg1_avg,1), size(cumLA_cfg1_avg,2), 3, nY);
        LA_cfg_hsv = zeros(size(LA_c_cfg1_avg,1), size(LA_c_cfg1_avg,2), 3, nY);
        PRRc = zeros(size(PhR_c_cfg1_avg,1), size(PhR_c_cfg1_avg,2), 3, nY, 'uint8');
        cumLA_Ms_cfg1_rmBG_hsv = zeros(size(cumLA_Ms_cfg1_rmBG,1), size(cumLA_Ms_cfg1_rmBG,2), 3, nY);
        LA_Ms_cfg1_rmBG_hsv = zeros(size(LA_Ms_cfg1_rmBG,1), size(LA_Ms_cfg1_rmBG,2), 3, nY);
        PRRc_rmBG = zeros(size(PhR_Ms_cfg1_rmBG,1), size(PhR_Ms_cfg1_rmBG,2), 3, nY, 'uint8');
        
        % 生成所有彩色编码
        for iY = 1:nY
            cumLA_cfg_hsv(:, :, :, iY) = quColoring(cumLA_cfg1_avg(:, :, :, iY), rotAngle);
            LA_cfg_hsv(:, :, :, iY) = quColoring(LA_c_cfg1_avg(:, :, :, iY), rotAngle);
            PRRc(:, :, :, iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_avg(:, :, iY), PRRrg) * 256), parula(256)) * 256);
            cumLA_Ms_cfg1_rmBG_hsv(:, :, :, iY) = quColoring(cumLA_Ms_cfg1_rmBG(:, :, :, iY), rotAngle);
            LA_Ms_cfg1_rmBG_hsv(:, :, :, iY) = quColoring(LA_Ms_cfg1_rmBG(:, :, :, iY), rotAngle);
            PRRc_rmBG(:, :, :, iY) = uint8(ind2rgb(uint8(mat2gray(PhR_Ms_cfg1_rmBG(:, :, iY), PRRrg) * 256), parula(256)) * 256);
        end
        
        % 应用边界mask：将边界以上的区域设为纯黑色（只处理前nZ_save层）
        for iY = 1:nY
            for iX = 1:size(cumLA_cfg_hsv, 2)
                surface_pos = round(topLines(iX, iY));
                
                % 确定需要涂黑的区域截止点
                mask_end_idx = 0;
                if surface_pos > nZ_save
                    mask_end_idx = nZ_save; % 涂黑整个深度
                elseif surface_pos > 1
                    mask_end_idx = surface_pos - 1; % 涂黑边界以上
                end
                
                if mask_end_idx > 0
                    % 对累积光轴HSV图像应用mask
                    cumLA_cfg_hsv(1:mask_end_idx, iX, :, iY) = 0;
                    % 对局部光轴HSV图像应用mask
                    LA_cfg_hsv(1:mask_end_idx, iX, :, iY) = 0;
                    % 对延迟相位彩色图像应用mask
                    PRRc(1:mask_end_idx, iX, :, iY) = 0;
                    % 对去背景后的累积光轴HSV图像应用mask
                    cumLA_Ms_cfg1_rmBG_hsv(1:mask_end_idx, iX, :, iY) = 0;
                    % 对去背景后的局部光轴HSV图像应用mask
                    LA_Ms_cfg1_rmBG_hsv(1:mask_end_idx, iX, :, iY) = 0;
                    % 对去背景后的延迟相位彩色图像应用mask
                    PRRc_rmBG(1:mask_end_idx, iX, :, iY) = 0;
                end
            end
        end
        
        % 写入 DCM 文件到 dcm 文件夹
        dicomwrite_with_lowquality(uint8(255 * (cumLA_cfg1_avg / 2 + 0.5)), fullfile(dcm_dir, [name, '_2-1_cumLA-cfg1-', num2str(nr), 'repAvg.dcm']), dcm_low_quality_dir, low_quality_scale);
        dicomwrite_with_lowquality(uint8(255 * cumLA_cfg_hsv), fullfile(dcm_dir, [name, '_2-2_cumLA-cfg1-', num2str(nr), 'repAvg_hsvColoring.dcm']), dcm_low_quality_dir, low_quality_scale);
        % dicomwrite(uint8(255 * (LA_c_cfg1_avg / 2 + 0.5)), fullfile(dcm_dir, [name, '_2-3_drLA-cfg1-', num2str(nr), 'repAvg.dcm']));
        % dicomwrite(uint8(255 * LA_cfg_hsv), fullfile(dcm_dir, [name, '_2-4_drLA-cfg1-', num2str(nr), 'repAvg_hsvColoring.dcm']));
        dicomwrite_with_lowquality(PRRc, fullfile(dcm_dir, [name, '_2-5_PhR-cfg1-', num2str(nr), 'repAvg.dcm']), dcm_low_quality_dir, low_quality_scale);
        % dicomwrite(uint8(255 * (cumLA_Ms_cfg1_rmBG / 2 + 0.5)), fullfile(dcm_dir, [name, '_2-6_cumLA_rmBG-cfg1-', num2str(nr), 'repAvg.dcm']));
        % dicomwrite(uint8(255 * cumLA_Ms_cfg1_rmBG_hsv), fullfile(dcm_dir, [name, '_2-7_cumLA_rmBG-cfg1-', num2str(nr), 'repAvg_hsvColoring.dcm']));
        % dicomwrite(uint8(255 * (LA_Ms_cfg1_rmBG / 2 + 0.5)), fullfile(dcm_dir, [name, '_2-8_drLA_rmBG-cfg1-', num2str(nr), 'repAvg.dcm']));
        % dicomwrite(uint8(255 * LA_Ms_cfg1_rmBG_hsv), fullfile(dcm_dir, [name, '_2-9_drLA_rmBG-cfg1-', num2str(nr), 'repAvg_hsvColoring.dcm']));
        % dicomwrite(PRRc_rmBG, fullfile(dcm_dir, [name, '_2-10_PhR_rmBG-cfg1-', num2str(nr), 'repAvg.dcm']));
        
        fprintf('彩色编码图像生成完成\n');
        % ensure tiff_dir exists (used later even if flatten disabled)
        tiff_dir = fullfile(foutputdir, 'tiff');
        if ~exist(tiff_dir, 'dir'), mkdir(tiff_dir); end
        
        % Print the current flag to help debugging
        fprintf('展平并生成 En-face 开关: %d\n', double(isfield(params.processing,'enable_flatten_enface') && params.processing.enable_flatten_enface));
        if ~isfield(params.processing,'enable_flatten_enface') || params.processing.enable_flatten_enface
            fprintf('\n========== 开始展平3D数据 ==========\n');
            fprintf('\n启用展平并生成 En-face (enable_flatten_enface = 1)\n');
            % tiff_dir 已在外部创建（以便在禁用展平时仍可执行 DCM->TIFF 转换）
        
        % 获取原始数据维度
        Struc_3D_orig = squeeze(Struc);  % [Z, X, Y]
        [nZ_orig, nX, nY] = size(Struc_3D_orig);
        [nZ_cumLA, ~, ~, ~] = size(cumLA_cfg_hsv);
        [nZ_PhR, ~, ~, ~] = size(PRRc);
        Avnum = params.polarization.Avnum;
        
        % (1) 展平 Struc - 基于 topLines 对齐表面
        fprintf('展平 Struc 数据...\n');
        t_flatten = tic;
        % 只使用前300层进行展平（减少20层）
        nZ_flatten = min(300, nZ_orig);
        Struc_flat = zeros(nZ_flatten, nX, nY, 'single');
        topLines_round = round(topLines);  % [nX, nY]
        
        for iy = 1:nY
            if mod(iy, max(1, floor(nY/5))) == 0 || iy == nY
                print_progress(iy, nY, 'Struc', 20);
            end
            for ix = 1:nX
                surf_z = topLines_round(ix, iy);
                
                % 如果表面位置超出范围，跳过展平
                if surf_z > nZ_orig
                    continue; 
                end
                
                surf_z = max(1, min(nZ_orig, surf_z));
                
                % 将表面对齐到固定深度（第1层）
                target_surf_z = 1;
                shift = target_surf_z - surf_z;
                
                % 复制并平移数据（只处理前300层）
                for z = 1:nZ_flatten
                    z_src = z - shift;
                    if z_src >= 1 && z_src <= nZ_orig
                        Struc_flat(z, ix, iy) = Struc_3D_orig(z_src, ix, iy);
                    end
                end
            end
        end
        fprintf('\n  Struc 展平完成，耗时: %.2f 秒\n', toc(t_flatten));
        
        % (2) 展平 cumLA_cfg_hsv
        fprintf('展平 cumLA 数据...\n');
        t_flatten = tic;
        % 对应Struc的300层，cumLA也使用相应层数
        nZ_cumLA_flatten = min(nZ_cumLA, nZ_flatten);
        cumLA_flat = zeros(nZ_cumLA_flatten, nX, 3, nY, 'single');
        
        for iy = 1:nY
            if mod(iy, max(1, floor(nY/5))) == 0 || iy == nY
                print_progress(iy, nY, 'cumLA', 20);
            end
            for ix = 1:nX
                surf_z = topLines_round(ix, iy);
                surf_z = max(1, min(nZ_orig, surf_z));
                
                target_surf_z = 1;
                shift = target_surf_z - surf_z;
                
                for z = 1:nZ_cumLA_flatten
                    z_src = z - shift;  % 修复：与Struc对齐，从顶端开始
                    if z_src >= 1 && z_src <= nZ_cumLA
                        cumLA_flat(z, ix, :, iy) = cumLA_cfg_hsv(z_src, ix, :, iy);
                    end
                end
            end
        end
        fprintf('\n  cumLA 展平完成，耗时: %.2f 秒\n', toc(t_flatten));
        
        % (3) 展平 PRRc (PhR)
        fprintf('展平 PhR 数据...\n');
        t_flatten = tic;
        % 对应Struc的300层，PhR也使用相应层数
        nZ_PhR_flatten = min(nZ_PhR, nZ_flatten);
        PhR_flat = zeros(nZ_PhR_flatten, nX, 3, nY, 'uint8');
        
        for iy = 1:nY
            if mod(iy, max(1, floor(nY/5))) == 0 || iy == nY
                print_progress(iy, nY, 'PhR', 20);
            end
            for ix = 1:nX
                surf_z = topLines_round(ix, iy);
                surf_z = max(1, min(nZ_orig, surf_z));
                
                target_surf_z = 1;
                shift = target_surf_z - surf_z;
                
                for z = 1:nZ_PhR_flatten
                    z_src = z - shift;  % 修复：与Struc对齐，从顶端开始
                    if z_src >= 1 && z_src <= nZ_PhR
                        PhR_flat(z, ix, :, iy) = PRRc(z_src, ix, :, iy);
                    end
                end
            end
        end
        fprintf('\n  PhR 展平完成，耗时: %.2f 秒\n', toc(t_flatten));
        
        %% ========== 步骤2a: 基于后巩膜边界展平数据（如果有） ==========
        if has_sclera_boundary
            fprintf('\n========== 基于后巩膜边界展平数据 ==========\n');
            
            % 使用后巩膜边界代替表面边界进行展平
            sclera_boundary_round = round(sclera_boundary_data);  % [nX, nY]
            
            % (1) 展平 Struc - 基于后巩膜边界对齐
            fprintf('基于后巩膜边界展平 Struc 数据...\n');
            t_sclera_flatten = tic;
            Struc_sclera_flat = zeros(nZ_flatten, nX, nY, 'single');
            
            for iy = 1:nY
                if mod(iy, max(1, floor(nY/5))) == 0 || iy == nY
                    print_progress(iy, nY, 'Struc_sclera', 20);
                end
                for ix = 1:nX
                    sclera_z = sclera_boundary_round(ix, iy);
                    sclera_z = max(1, min(nZ_orig, sclera_z));
                    
                    % 将后巩膜边界对齐到固定深度（第1层）
                    target_sclera_z = 1;
                    shift = target_sclera_z - sclera_z;
                    
                    % 复制并平移数据
                    for z = 1:nZ_flatten
                        z_src = z - shift;
                        if z_src >= 1 && z_src <= nZ_orig
                            Struc_sclera_flat(z, ix, iy) = Struc_3D_orig(z_src, ix, iy);
                        end
                    end
                end
            end
            fprintf('\n  Struc (后巩膜) 展平完成，耗时: %.2f 秒\n', toc(t_sclera_flatten));
            
            % (2) 展平 cumLA - 基于后巩膜边界
            fprintf('基于后巩膜边界展平 cumLA 数据...\n');
            t_sclera_cumla = tic;
            cumLA_sclera_flat = zeros(nZ_cumLA_flatten, nX, 3, nY, 'single');
            
            for iy = 1:nY
                if mod(iy, max(1, floor(nY/5))) == 0 || iy == nY
                    print_progress(iy, nY, 'cumLA_sclera', 20);
                end
                for ix = 1:nX
                    sclera_z = sclera_boundary_round(ix, iy);
                    sclera_z = max(1, min(nZ_orig, sclera_z));
                    
                    target_sclera_z = 1;
                    shift = target_sclera_z - sclera_z;
                    
                    for z = 1:nZ_cumLA_flatten
                        z_src = z - shift;  % 修复：与Struc对齐，从顶端开始
                        if z_src >= 1 && z_src <= nZ_cumLA
                            cumLA_sclera_flat(z, ix, :, iy) = cumLA_cfg_hsv(z_src, ix, :, iy);
                        end
                    end
                end
            end
            fprintf('\n  cumLA (后巩膜) 展平完成，耗时: %.2f 秒\n', toc(t_sclera_cumla));
            
            % (3) 展平 PhR - 基于后巩膜边界
            fprintf('基于后巩膜边界展平 PhR 数据...\n');
            t_sclera_phr = tic;
            PhR_sclera_flat = zeros(nZ_PhR_flatten, nX, 3, nY, 'uint8');
            
            for iy = 1:nY
                if mod(iy, max(1, floor(nY/5))) == 0 || iy == nY
                    print_progress(iy, nY, 'PhR_sclera', 20);
                end
                for ix = 1:nX
                    sclera_z = sclera_boundary_round(ix, iy);
                    sclera_z = max(1, min(nZ_orig, sclera_z));
                    
                    target_sclera_z = 1;
                    shift = target_sclera_z - sclera_z;
                    
                    for z = 1:nZ_PhR_flatten
                        z_src = z - shift;  % 修复：与Struc对齐，从顶端开始
                        if z_src >= 1 && z_src <= nZ_PhR
                            PhR_sclera_flat(z, ix, :, iy) = PRRc(z_src, ix, :, iy);
                        end
                    end
                end
            end
            fprintf('\n  PhR (后巩膜) 展平完成，耗时: %.2f 秒\n', toc(t_sclera_phr));
            fprintf('========== 后巩膜边界展平完成 ==========\n\n');
        end
        
        %% ========== 步骤2: 保存展平后的 DCM（生成新文件，命名中加上 flat） ==========
        fprintf('\n保存展平后的 DCM 文件...\n');
        
        % 保存展平后的 Struc（新文件，加上 flat 标识）
        Struc_flat_4d = permute(Struc_flat, [1, 2, 4, 3]);  % [Z, X, 1, Y] for DCM
        dicomwrite_with_lowquality(uint8(255 * Struc_flat_4d), fullfile(dcm_dir, [name, '_1-1_Struc_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        
        % 保存展平后的 cumLA（新文件，加上 flat 标识）
        dicomwrite_with_lowquality(uint8(255 * cumLA_flat), fullfile(dcm_dir, [name, '_2-2_cumLA-cfg1-', num2str(nr), 'repAvg_hsvColoring_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        
        % 保存展平后的 PhR（新文件，加上 flat 标识）
        dicomwrite_with_lowquality(PhR_flat, fullfile(dcm_dir, [name, '_2-5_PhR-cfg1-', num2str(nr), 'repAvg_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        
        % 如果有后巩膜边界展平数据，也保存这些文件
        if has_sclera_boundary
            fprintf('保存后巩膜边界展平的 DCM 文件...\n');
            
            % 保存后巩膜展平的 Struc
            Struc_sclera_flat_4d = permute(Struc_sclera_flat, [1, 2, 4, 3]);
            dicomwrite_with_lowquality(uint8(255 * Struc_sclera_flat_4d), fullfile(dcm_dir, [name, '_1-1_Struc_sclera_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
            
            % 保存后巩膜展平的 cumLA
            dicomwrite_with_lowquality(uint8(255 * cumLA_sclera_flat), fullfile(dcm_dir, [name, '_2-2_cumLA-cfg1-', num2str(nr), 'repAvg_hsvColoring_sclera_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
            
            % 保存后巩膜展平的 PhR
            dicomwrite_with_lowquality(PhR_sclera_flat, fullfile(dcm_dir, [name, '_2-5_PhR-cfg1-', num2str(nr), 'repAvg_sclera_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
            
            fprintf('后巩膜边界展平的 DCM 保存完成\n');
        end
        
        fprintf('展平后的 DCM 保存完成\n');
        
        % 释放不再需要的中间数组（保留cumLA_cfg_hsv和PRRc供后续En-face生成使用）
        fprintf('释放中间变量以降低内存占用...\n');
        clear LA_c_cfg1_avg PhR_c_cfg1_avg LA_c_cfg1_eig;
        clear LA_Ms_cfg1_rmBG PhR_Ms_cfg1_rmBG cumLA_Ms_cfg1_rmBG;
        clear Smap_avg Smap_rep1;
    
        %% ========== 步骤3: 从展平后的数据直接切片生成 En-face ==========
        
        % 使用展平后的数据维度
        [nZ_out, nX_out, nY_out] = size(Struc_flat);

    % % 1. Orthogonal B-scan (Y-Z plane, 沿着慢轴和深度轴)
    % % 选取 X 轴(快轴)中间位置切片
    % x_mid = round(nX_out / 2);
    % if x_mid < 1, x_mid = 1; end
    % ortho_bscan = squeeze(Struc_3D(:, x_mid, :)); % [Z, Y]
    % dicomwrite(uint8(255 * ortho_bscan), fullfile(foutputdir, [name, '_3-1_Ortho_Bscan_X', num2str(x_mid), '.dcm']));

    % % 2. En-face (C-scan, X-Y plane, 平行于样品表面)
    % % 生成所有深度的结构图像 En-face 切面，多帧 DCM
    % fprintf('生成结构图像 En-face 多帧 DCM...\n');
    % en_face_stack = zeros(nY_out, nX_out, 1, nZ_out, 'uint8'); % [Y, X, 1, Z] for multi-frame
    % for z = 1:nZ_out
    %     % 提取当前深度的切片 [X, Y]，然后转置为 [Y, X]
    %     en_face = squeeze(Struc_3D(z, :, :))'; % [Y, X]
    %     en_face_stack(:, :, 1, z) = uint8(255 * en_face);
    % end
    % dicomwrite(en_face_stack, fullfile(foutputdir, [name, '_3-2_Enface_Struc_MultiFrame.dcm']));

    % (c) 从展平数据直接切片生成 En-face
    if exist('topLines', 'var') && ~isempty(topLines)
            fprintf('\n========== 开始生成 En-face 切片 ==========\n');
            fprintf('从展平数据生成 En-face 切片...\n');
            t_enface = tic;
            
            % 预分配 DCM 数组
            enface_struc_stack = zeros(nY_out, nX_out, 1, nZ_out, 'uint8');
            
            % 直接从展平数据切片（无需再计算偏移）
            for z = 1:nZ_out
                if mod(z, max(1, floor(nZ_out/5))) == 0 || z == nZ_out
                    print_progress(z, nZ_out, 'En-face', 20);
                end
                
                % 从展平数据直接提取 [X, Y] 切片，转置为 [Y, X]
                en_face_slice = squeeze(Struc_flat(z, :, :))';  % [Y, X]
                
                % 保存到 DCM 数组
                enface_struc_stack(:, :, 1, z) = uint8(255 * en_face_slice);
            end
            
            % 保存 En-face DCM
            dicomwrite_with_lowquality(enface_struc_stack, fullfile(dcm_dir, [name, '_3-3_Enface_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
            fprintf('\n  Struc En-face 耗时: %.2f 秒\n', toc(t_enface));
            
            % 释放临时变量（保留Struc等原始数据供后续使用）
            clear enface_struc_stack en_face_slice Struc_3D_orig Struc_flat;
    end

    % % 3. DOPU En-face (更可能显示血管)
    % fprintf('生成 DOPU En-face 多帧 DCM...\n');
    % dopu_3D = dopu_splitSpectrum; % [Z, X, Y]
    % dopu_enface_stack = zeros(nY_out, nX_out, 1, nZ_out, 'uint8');
    % for z = 1:nZ_out
    %     dopu_enface = squeeze(dopu_3D(z, :, :))'; % [Y, X]
    %     dopu_enface_stack(:, :, 1, z) = uint8(255 * dopu_enface);
    % end
    % dicomwrite(dopu_enface_stack, fullfile(foutputdir, [name, '_3-4_Enface_DOPU_MultiFrame.dcm']));
    % dopu_with_boundary = dopu_splitSpectrum;
    % for iY = 1:size(dopu_with_boundary, 3)
    %     for iX = 1:size(dopu_with_boundary, 2)
    %         surface_pos = round(topLines(iX, iY));
    %         if surface_pos > 1 && surface_pos <= size(dopu_with_boundary, 1)
    %             dopu_with_boundary(1:surface_pos-1, iX, iY) = 0;  % 边界以上的部分设为黑色
    %         end
    %     end
    % end

    % 4. cumLA En-face - 从展平数据直接切片
    if exist('topLines', 'var') && ~isempty(topLines)
        fprintf('从展平 cumLA 数据生成 En-face...\n');
        t_cumla = tic;
        [nZ_cumLA, nX_cumLA, ~, nY_cumLA] = size(cumLA_flat);

        % 预分配 DCM 数组
        enface_cumla_stack = zeros(nY_cumLA, nX_cumLA, 3, nZ_cumLA, 'uint8');
        
        % 直接从展平数据切片
        for z = 1:nZ_cumLA
            if mod(z, max(1, floor(nZ_cumLA/5))) == 0 || z == nZ_cumLA
                print_progress(z, nZ_cumLA, 'cumLA', 20);
            end
            
            % 从展平数据直接提取 [X, 3, Y] 切片，重排为 [Y, X, 3]
            slice_rgb = permute(squeeze(cumLA_flat(z, :, :, :)), [3, 1, 2]);  % [Y, X, 3]
            
            % 保存到 DCM 数组
            enface_cumla_stack(:, :, :, z) = uint8(255 * slice_rgb);
        end
        
        % 保存 En-face DCM
        dicomwrite_with_lowquality(enface_cumla_stack, fullfile(dcm_dir, [name, '_3-5_Enface_cumLA_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  cumLA En-face 耗时: %.2f 秒\n', toc(t_cumla));
        
        % 释放大数组
        clear cumLA_flat enface_cumla_stack;
    end

    % 5. PhR En-face - 从展平数据直接切片
    if exist('topLines', 'var') && ~isempty(topLines)
        fprintf('从展平 PhR 数据生成 En-face...\n');
        t_phr = tic;
        [nZ_PhR, nX_PhR, ~, nY_PhR] = size(PhR_flat);

        % 预分配 DCM 数组
        enface_phr_stack = zeros(nY_PhR, nX_PhR, 3, nZ_PhR, 'uint8');
        
        % 直接从展平数据切片
        for z = 1:nZ_PhR
            if mod(z, max(1, floor(nZ_PhR/5))) == 0 || z == nZ_PhR
                print_progress(z, nZ_PhR, 'PhR', 20);
            end
            
            % 从展平数据直接提取 [X, 3, Y] 切片，重排为 [Y, X, 3]
            slice_rgb = permute(squeeze(PhR_flat(z, :, :, :)), [3, 1, 2]);  % [Y, X, 3]
            
            % 保存到 DCM 数组
            enface_phr_stack(:, :, :, z) = slice_rgb;
        end
        
        % 保存 En-face DCM
        dicomwrite_with_lowquality(enface_phr_stack, fullfile(dcm_dir, [name, '_3-6_Enface_PhR_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  PhR En-face 耗时: %.2f 秒\n', toc(t_phr));
        
        % 释放大数组
        clear PhR_flat PRRc_rmBG LA_cfg_hsv LA_Ms_cfg1_rmBG_hsv enface_phr_stack;
    end
    
    %% ========== 步骤4: 从后巩膜边界展平数据生成 En-face（如果有）==========
    if has_sclera_boundary
        fprintf('\n========== 基于后巩膜边界生成 En-face 切片 ==========\n');
        
        % 6. Struc En-face - 基于后巩膜边界展平
        fprintf('从后巩膜展平 Struc 数据生成 En-face...\n');
        t_struc_sclera = tic;
        [nZ_sclera_struc, nX_sclera_struc, nY_sclera_struc] = size(Struc_sclera_flat);
        
        enface_struc_sclera_stack = zeros(nY_sclera_struc, nX_sclera_struc, 1, nZ_sclera_struc, 'uint8');
        
        for z = 1:nZ_sclera_struc
            if mod(z, max(1, floor(nZ_sclera_struc/5))) == 0 || z == nZ_sclera_struc
                print_progress(z, nZ_sclera_struc, 'Struc_sclera', 20);
            end
            
            % 从展平数据直接提取 [X, Y] 切片，转置为 [Y, X]
            en_face_slice = squeeze(Struc_sclera_flat(z, :, :))';  % [Y, X]
            enface_struc_sclera_stack(:, :, 1, z) = uint8(255 * en_face_slice);
        end
        
        dicomwrite_with_lowquality(enface_struc_sclera_stack, fullfile(dcm_dir, [name, '_4-1_Enface_Struc_sclera_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  Struc (后巩膜) En-face 耗时: %.2f 秒\n', toc(t_struc_sclera));
        clear enface_struc_sclera_stack Struc_sclera_flat;
        
        % 7. cumLA En-face - 基于后巩膜边界展平
        fprintf('从后巩膜展平 cumLA 数据生成 En-face...\n');
        t_cumla_sclera = tic;
        [nZ_sclera_cumLA, nX_sclera_cumLA, ~, nY_sclera_cumLA] = size(cumLA_sclera_flat);
        
        enface_cumla_sclera_stack = zeros(nY_sclera_cumLA, nX_sclera_cumLA, 3, nZ_sclera_cumLA, 'uint8');
        
        for z = 1:nZ_sclera_cumLA
            if mod(z, max(1, floor(nZ_sclera_cumLA/5))) == 0 || z == nZ_sclera_cumLA
                print_progress(z, nZ_sclera_cumLA, 'cumLA_sclera', 20);
            end
            
            % 从展平数据直接提取 [X, 3, Y] 切片，重排为 [Y, X, 3]
            slice_rgb = permute(squeeze(cumLA_sclera_flat(z, :, :, :)), [3, 1, 2]);  % [Y, X, 3]
            enface_cumla_sclera_stack(:, :, :, z) = uint8(255 * slice_rgb);
        end
        
        dicomwrite_with_lowquality(enface_cumla_sclera_stack, fullfile(dcm_dir, [name, '_4-2_Enface_cumLA_sclera_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  cumLA (后巩膜) En-face 耗时: %.2f 秒\n', toc(t_cumla_sclera));
        clear enface_cumla_sclera_stack cumLA_sclera_flat;
        
        % 8. PhR En-face - 基于后巩膜边界展平
        fprintf('从后巩膜展平 PhR 数据生成 En-face...\n');
        t_phr_sclera = tic;
        [nZ_sclera_PhR, nX_sclera_PhR, ~, nY_sclera_PhR] = size(PhR_sclera_flat);
        
        enface_phr_sclera_stack = zeros(nY_sclera_PhR, nX_sclera_PhR, 3, nZ_sclera_PhR, 'uint8');
        
        for z = 1:nZ_sclera_PhR
            if mod(z, max(1, floor(nZ_sclera_PhR/5))) == 0 || z == nZ_sclera_PhR
                print_progress(z, nZ_sclera_PhR, 'PhR_sclera', 20);
            end
            
            % 从展平数据直接提取 [X, 3, Y] 切片，重排为 [Y, X, 3]
            slice_rgb = permute(squeeze(PhR_sclera_flat(z, :, :, :)), [3, 1, 2]);  % [Y, X, 3]
            enface_phr_sclera_stack(:, :, :, z) = slice_rgb;
        end
        
        dicomwrite_with_lowquality(enface_phr_sclera_stack, fullfile(dcm_dir, [name, '_4-3_Enface_PhR_sclera_flat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  PhR (后巩膜) En-face 耗时: %.2f 秒\n', toc(t_phr_sclera));
        clear enface_phr_sclera_stack PhR_sclera_flat;
        
        fprintf('========== 后巩膜边界 En-face 生成完成 ==========\n');
    else
        fprintf('\n跳过展平及 En-face 生成功能 (enable_flatten_enface = 0)\n');
    end % end of enable_flatten_enface check

    %% ========== 生成非展平En-face切片（直接从原始数据切片） ==========
    if isfield(params.processing,'enable_enface_noflat') && params.processing.enable_enface_noflat
        fprintf('\n========== 开始生成非展平En-face切片 ==========\n');
        fprintf('启用非展平En-face生成 (enable_enface_noflat = 1)\n');
        
        % 使用前300层原始数据生成En-face
        Struc_noflat = Struc(1:nZ_save, :, :, :);  % [Z, X, 1, Y]
        [nZ_noflat, nX_noflat, ~, nY_noflat] = size(Struc_noflat);
        
        % 1. Struc En-face - 非展平版本
        fprintf('从原始 Struc 数据生成 En-face...\n');
        t_enface_noflat = tic;
        enface_struc_noflat_stack = zeros(nY_noflat, nX_noflat, 1, nZ_noflat, 'uint8');
        
        for z = 1:nZ_noflat
            if mod(z, max(1, floor(nZ_noflat/5))) == 0 || z == nZ_noflat
                print_progress(z, nZ_noflat, 'Struc', 20);
            end
            % 直接从原始数据提取 [X, Y] 切片，转置为 [Y, X]
            en_face_slice = squeeze(Struc_noflat(z, :, 1, :))';  % [Y, X]
            enface_struc_noflat_stack(:, :, 1, z) = uint8(255 * en_face_slice);
        end
        
        dicomwrite_with_lowquality(enface_struc_noflat_stack, fullfile(dcm_dir, [name, '_3-3_Enface_noflat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  Struc En-face (非展平) 耗时: %.2f 秒\n', toc(t_enface_noflat));
        clear enface_struc_noflat_stack;
        
        % 2. cumLA En-face - 非展平版本
        fprintf('从原始 cumLA 数据生成 En-face...\n');
        t_cumla_noflat = tic;
        
        % cumLA_cfg_hsv 已经是彩色编码后的数据 [Z, X, 3, Y]
        [nZ_cumLA_orig, nX_cumLA_orig, ~, nY_cumLA_orig] = size(cumLA_cfg_hsv);
        nZ_cumLA_save = min(nZ_save, nZ_cumLA_orig);  % 使用前300层
        
        enface_cumla_noflat_stack = zeros(nY_cumLA_orig, nX_cumLA_orig, 3, nZ_cumLA_save, 'uint8');
        
        for z = 1:nZ_cumLA_save
            if mod(z, max(1, floor(nZ_cumLA_save/5))) == 0 || z == nZ_cumLA_save
                print_progress(z, nZ_cumLA_save, 'cumLA', 20);
            end
            % 从原始数据直接提取 [X, 3, Y] 切片，重排为 [Y, X, 3]
            slice_rgb = permute(squeeze(cumLA_cfg_hsv(z, :, :, :)), [3, 1, 2]);  % [Y, X, 3]
            enface_cumla_noflat_stack(:, :, :, z) = uint8(255 * slice_rgb);
        end
        
        dicomwrite_with_lowquality(enface_cumla_noflat_stack, fullfile(dcm_dir, [name, '_3-5_Enface_cumLA_noflat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  cumLA En-face (非展平) 耗时: %.2f 秒\n', toc(t_cumla_noflat));
        clear enface_cumla_noflat_stack;
        
        % 3. PhR En-face - 非展平版本
        fprintf('从原始 PhR 数据生成 En-face...\n');
        t_phr_noflat = tic;
        
        % PRRc 已经是彩色编码后的数据 [Z, X, 3, Y]
        [nZ_PhR_orig, nX_PhR_orig, ~, nY_PhR_orig] = size(PRRc);
        nZ_PhR_save = min(nZ_save, nZ_PhR_orig);  % 使用前300层
        
        enface_phr_noflat_stack = zeros(nY_PhR_orig, nX_PhR_orig, 3, nZ_PhR_save, 'uint8');
        
        for z = 1:nZ_PhR_save
            if mod(z, max(1, floor(nZ_PhR_save/5))) == 0 || z == nZ_PhR_save
                print_progress(z, nZ_PhR_save, 'PhR', 20);
            end
            % 从原始数据直接提取 [X, 3, Y] 切片，重排为 [Y, X, 3]
            slice_rgb = permute(squeeze(PRRc(z, :, :, :)), [3, 1, 2]);  % [Y, X, 3]
            enface_phr_noflat_stack(:, :, :, z) = slice_rgb;
        end
        
        dicomwrite_with_lowquality(enface_phr_noflat_stack, fullfile(dcm_dir, [name, '_3-6_Enface_PhR_noflat.dcm']), dcm_low_quality_dir, low_quality_scale);
        fprintf('\n  PhR En-face (非展平) 耗时: %.2f 秒\n', toc(t_phr_noflat));
        clear enface_phr_noflat_stack Struc_noflat;
        
        fprintf('\n非展平En-face切片生成完成\n');
    else
        fprintf('\n跳过非展平En-face生成 (enable_enface_noflat = 0)\n');
    end % end of enable_enface_noflat check
    
    %% 清理所有En-face相关的原始数据
    clear Struc Struc_with_boundary cumLA_cfg_hsv PRRc;

    end % save results

    %% DCM到TIFF转换功能
    if params.tiff.saveDicom && params.tiff.make_tiff
        fprintf('\n开始DCM到TIFF转换...\n');
        try
            % 调用本地DCM到TIFF转换函数，从 dcm 子文件夹读取，输出到 tiff 子文件夹
            convert_dcm_to_tiff_local(dcm_dir, params.tiff.tiff_frame, tiff_dir);
        catch ME
            fprintf('DCM到TIFF转换出错: %s\n', ME.message);
            % 转换出错不影响主流程，继续执行
        end
    end

    fprintf('文件 %s 处理完成!\n', display_name);
    fprintf('输出目录: %s\n', foutputdir);
    
    % 【内存监控】显示最终内存状态
    try
        if ispc
            m = memory;
            fprintf('最终可用内存: %.1f GB\n', m.MemAvailableAllArrays/1024^3);
        end
    catch
        % 忽略内存监控错误
    end

    % 显示处理时间统计
    file_proc_time = toc(file_start_time);
    fprintf('处理时间: %.2f 秒 (%.2f 分钟)\n', file_proc_time, file_proc_time/60);
    fprintf('\n');

    % 清理并行池（在函数结束时）
    % 仅在本函数启动了并行池且配置要求自动关闭时才删除（默认保持并行池存活以提高效率）
    if exist('poolStartedHere','var') && poolStartedHere && isfield(params,'parallel') && isfield(params.parallel,'autoClosePool') && params.parallel.autoClosePool
        fprintf('关闭并行池（由本函数启动且 autoClosePool = true）\n');
        delete(gcp('nocreate'));
    else
        fprintf('保留并行池（poolStartedHere=%d, autoClosePool=%d）\n', double(exist('poolStartedHere','var') && poolStartedHere), double(isfield(params,'parallel') && isfield(params.parallel,'autoClosePool') && params.parallel.autoClosePool));
    end
    end
end
return;
