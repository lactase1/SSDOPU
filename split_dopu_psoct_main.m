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
data_path   = 'C:\yongxin.wang/Data/Process_Data/Optic_Disc_rep';
% σ * 6 + 1 // σ * 4 + 1
output_base = 'C:\yongxin.wang/Data/Process_Data/Optic_Disc_rep/Output-dopu-adj/dopu_13layer_3_13';
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
    
    % 逐个处理文件
    for i = 1:length(oct_files)
        full_path = fullfile(data_path, oct_files(i).name);
        fprintf('\n==================================================\n');
        fprintf('开始处理第 %d/%d 个文件: %s\n', i, length(oct_files), oct_files(i).name);
        fprintf('==================================================\n');
        
        try
            % 记录开始时间
            tic;
            
            % 调用单文件处理函数
            rPSOCT_process_single_file(full_path, output_base);
            
            % 计算并显示处理时间
            proc_time = toc;
            fprintf('文件 %s 处理成功, 耗时: %.2f 秒 (%.2f 分钟)\n', oct_files(i).name, proc_time, proc_time/60);
        catch ME
            fprintf('处理文件 %s 时出错: %s\n', oct_files(i).name, ME.message);
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
        [~, outbase_name, ~] = fileparts(output_base);
        updates = struct();
        wrote_any = false;
        
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
        foutputdir = output_base;
        if ~exist(foutputdir, 'dir'), mkdir(foutputdir); end

        czrg = 1:320; % set z range
        topLines = ones(nX, nY);
        if params.processing.hasSeg
            matFilePath = fullfile(filepath, [name, '.mat']);
            if exist(matFilePath, 'file')
                try
                    load(matFilePath, 'topLines');
                    topLines = double(topLines);
                    topLines(topLines <= 1) = 1;
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
                        local_enableOutputAdaptive, local_kRL_output, local_kRU_output);
                
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
            
            % ========== 步骤4: 清空批次数据（释放内存）==========
            clear batch_Bs1 batch_Bs2 batch_dopu batch_Strus batch_Smap_avg batch_Smap_rep1;
            clear batch_topLines batch_LA_c_cfg1_avg batch_PhR_c_cfg1_avg batch_cumLA_cfg1_avg;
            clear batch_LA_Ms_cfg1_rmBG batch_PhR_Ms_cfg1_rmBG batch_cumLA_Ms_cfg1_rmBG;
            
        end % 分批循环结束
        
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
        dicomwrite(uint8(255 * (Struc)), fullfile(dcm_dir, [name, '_1-1_Struc.dcm']));

        % 创建带边界的结构图像副本（只保留前300层）
        Struc_with_boundary = Struc(1:nZ_save, :, :, :);
        % 对Struc应用边界处理：边界以上的部分变成黑色
        for iY = 1:size(Struc_with_boundary, 4)
            for iX = 1:size(Struc_with_boundary, 2)
                surface_pos = round(topLines(iX, iY));
                if surface_pos > 1 && surface_pos <= nZ_save
                    Struc_with_boundary(1:surface_pos-1, iX, :, iY) = 0;  % 边界以上的部分设为黑色
                end
            end
        end
        dicomwrite(uint8(255 * (Struc(1:nZ_save, :, :, :))), fullfile(dcm_dir, [name, '_1-1_Struc.dcm']));
        dicomwrite(uint8(255 * (Struc_with_boundary)), fullfile(dcm_dir, [name, '_1-2_Struc_with_boundary.dcm']));
        % dicomwrite(uint8(255 * (Smap_rep1 / 2 + 0.5)), fullfile(dcm_dir, [name, '_1-3_1rep-Stokes.dcm']));
        dicomwrite(uint8(255 * (Smap_avg(1:nZ_save, :, :, :) / 2 + 0.5)), fullfile(dcm_dir, [name, '_1-3_4rep-Stokes.dcm']));

        % 对1-4 DOPU图像应用阈值过滤（只保留前300层）
        dopu_thresholded = dopu_splitSpectrum(1:nZ_save, :, :);
        dopu_thresholded(dopu_thresholded <= 0.38) = 0;  % 小于等于0.5的设为0

        dicomwrite(uint8(255 * (permute(dopu_thresholded, [1 2 4 3]))), fullfile(dcm_dir, [name, '_1-4_dopu_SS.dcm']));
        
        if ~params.processing.hasSeg
            save(fullfile(filepath, [name, '.mat']), 'topLines', 'czrg');
        end
        rotAngle = 440;
        % 仅输出 cfg1 avg 相关结果（当前仅支持 avg 路径）
        PRRrg = [0 0.5];
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
        
        % 写入 DCM 文件到 dcm 文件夹
        dicomwrite(uint8(255 * (cumLA_cfg1_avg / 2 + 0.5)), fullfile(dcm_dir, [name, '_2-1_cumLA-cfg1-', num2str(nr), 'repAvg.dcm']));
        dicomwrite(uint8(255 * cumLA_cfg_hsv), fullfile(dcm_dir, [name, '_2-2_cumLA-cfg1-', num2str(nr), 'repAvg_hsvColoring.dcm']));
        % dicomwrite(uint8(255 * (LA_c_cfg1_avg / 2 + 0.5)), fullfile(dcm_dir, [name, '_2-3_drLA-cfg1-', num2str(nr), 'repAvg.dcm']));
        % dicomwrite(uint8(255 * LA_cfg_hsv), fullfile(dcm_dir, [name, '_2-4_drLA-cfg1-', num2str(nr), 'repAvg_hsvColoring.dcm']));
        dicomwrite(PRRc, fullfile(dcm_dir, [name, '_2-5_PhR-cfg1-', num2str(nr), 'repAvg.dcm']));
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
                    z_src = z - shift - Avnum;  % cumLA 从表面下 Avnum 开始
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
                    z_src = z - shift - Avnum;
                    if z_src >= 1 && z_src <= nZ_PhR
                        PhR_flat(z, ix, :, iy) = PRRc(z_src, ix, :, iy);
                    end
                end
            end
        end
        fprintf('\n  PhR 展平完成，耗时: %.2f 秒\n', toc(t_flatten));
        
        %% ========== 步骤2: 保存展平后的 DCM（生成新文件，命名中加上 flat） ==========
        fprintf('\n保存展平后的 DCM 文件...\n');
        
        % 保存展平后的 Struc（新文件，加上 flat 标识）
        Struc_flat_4d = permute(Struc_flat, [1, 2, 4, 3]);  % [Z, X, 1, Y] for DCM
        dicomwrite(uint8(255 * Struc_flat_4d), fullfile(dcm_dir, [name, '_1-1_Struc_flat.dcm']));
        
        % 保存展平后的 cumLA（新文件，加上 flat 标识）
        dicomwrite(uint8(255 * cumLA_flat), fullfile(dcm_dir, [name, '_2-2_cumLA-cfg1-', num2str(nr), 'repAvg_hsvColoring_flat.dcm']));
        
        % 保存展平后的 PhR（新文件，加上 flat 标识）
        dicomwrite(PhR_flat, fullfile(dcm_dir, [name, '_2-5_PhR-cfg1-', num2str(nr), 'repAvg_flat.dcm']));
        
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
            dicomwrite(enface_struc_stack, fullfile(dcm_dir, [name, '_3-3_Enface_flat.dcm']));
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
        dicomwrite(enface_cumla_stack, fullfile(dcm_dir, [name, '_3-5_Enface_cumLA_flat.dcm']));
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
        dicomwrite(enface_phr_stack, fullfile(dcm_dir, [name, '_3-6_Enface_PhR_flat.dcm']));
        fprintf('\n  PhR En-face 耗时: %.2f 秒\n', toc(t_phr));
        
        % 释放大数组
        clear PhR_flat PRRc_rmBG LA_cfg_hsv LA_Ms_cfg1_rmBG_hsv enface_phr_stack;
    end
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
        
        dicomwrite(enface_struc_noflat_stack, fullfile(dcm_dir, [name, '_3-3_Enface_noflat.dcm']));
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
        
        dicomwrite(enface_cumla_noflat_stack, fullfile(dcm_dir, [name, '_3-5_Enface_cumLA_noflat.dcm']));
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
        
        dicomwrite(enface_phr_noflat_stack, fullfile(dcm_dir, [name, '_3-6_Enface_PhR_noflat.dcm']));
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
return;
