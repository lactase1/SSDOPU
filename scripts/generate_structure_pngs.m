% generate_structure_pngs.m
% 仅提取结构（Struc）并把每个帧（B-scan）的结构图保存为 PNG
% 使用方法: 直接运行此脚本或在命令行中修改 data_path/output_base 等变量后运行

% 添加 function 目录到路径
script_dir = fileparts(mfilename('fullpath'));
function_path = fullfile(script_dir, '..', 'function');
if exist(function_path, 'dir')
    addpath(function_path);
end

%% 用户配置
data_path = 'C:\yongxin.wang\Data\Select_Data\Enface_Data\Jian'; % 修改为你的数据路径
output_base = 'D:\1-Liu Jian\yongxin.wang\Output\struct_pngs_Jian'; % 输出根目录
max_frames = -1; % <=0 表示全部帧
czrg_limit = 320; % 深度范围上限
verbose = true;

% 并行 / 性能相关配置
use_parallel = true;                 % 是否启用并行处理（parpool + parfor）
num_workers = min(48, feature('numcores')); % 可调整为你的可用 CPU 核心数（默认最多 48）
sample_frames = 200;                % 用于估算强度范围的采样帧数量（越小越快）
use_percentile = false;             % 是否使用分位数而非 min/max (如果 true 会采样更多帧并计算分位)
percentiles = [1, 99];             % 使用分位数时的百分位
use_gpu = false;                    % 是否尝试在每个 worker 上使用 GPU（须保证函数支持 gpuArray）

% 启动并行池（如果需要）
if use_parallel
    try
        if isempty(gcp('nocreate'))
            parpool('local', num_workers);
            fprintf('已启动并行池: %d workers\n', num_workers);
        else
            pool = gcp();
            fprintf('并行池已存在: %d workers\n', pool.NumWorkers);
        end
    catch ME
        warning('启动并行池失败: %s\n将退回为串行执行', ME.message);
        use_parallel = false;
    end
end

% 查找 oct 文件
oct_files = dir(fullfile(data_path, '*.oct'));
if isempty(oct_files)
    fprintf('未找到任何 .oct 文件: %s\n', data_path);
    return;
end

for fi = 1:length(oct_files)
    infile = fullfile(data_path, oct_files(fi).name);
    [~, name, ~] = fileparts(infile);
    fprintf('\n==== 处理 %s (%d/%d) ====', name, fi, length(oct_files));

    % 为该文件创建输出文件夹
    outdir = fullfile(output_base, name, 'struct_pngs');
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    fid = fopen(infile, 'r');
    if fid == -1
        warning('无法打开文件: %s', infile);
        continue;
    end

    % 读取头部（与主脚本一致）
    bob               = fread(fid, 1, 'uint32');
    SPL               = fread(fid, 1, 'double');
    nX                = fread(fid, 1, 'uint32'); %% number of Alines
    nY_raw            = fread(fid, 1, 'uint32'); %% number of B-scans
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

    nR = frame_per_pos;
    nY = floor(nY_raw / nR);
    if max_frames > 0
        nY = min(nY, max_frames);
    end

    % 深度范围
    czrg = 1:min(czrg_limit, floor(Blength/2));

    winG_whole = tukeywin(Blength, 0.25);

    % 参考背景
    Ref_ch1 = repmat(Bck1, [1, nX]);
    Ref_ch2 = repmat(Bck2, [1, nX]);

    % 第一遍：扫描所有帧以获取强度范围（min/max）用于一致归一化
    fprintf('\n第一遍扫描 (计算 min/max)...\n');
    global_min = inf;
    global_max = -inf;
    for iy = 1:nY
        if mod(iy, max(1, floor(nY/10))) == 0 || iy == nY
            fprintf('  扫描帧 %d/%d\n', iy, nY);
        end
        % 读取该帧的数据
        fseek(fid, bob + (SPL * nX * 2 + 2) * 2 * (nR * (iy - 1)), 'bof');
        Bs1 = zeros(SPL, nX, nR, 'single');
        Bs2 = zeros(SPL, nX, nR, 'single');
        for Ic = 1:nR
            fseek(fid, 4, 'cof');
            B = single(reshape(fread(fid, SPL * nX * 2, 'int16'), [SPL, nX * 2]));
            Bs1(:, :, Ic) = (B(:, 1:nX) - single(Ref_ch1));
            Bs2(:, :, Ic) = (B(:, nX + 1:end) - single(Ref_ch2));
        end

        % 计算 whole spectrum 图像并截取深度范围
        Bimg1_wholeStr = fft(Bs1 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
        Bimg2_wholeStr = fft(Bs2 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
        IMG1_wholeStr = Bimg1_wholeStr(czrg, :, :);
        IMG2_wholeStr = Bimg2_wholeStr(czrg, :, :);

        % 计算 Stokes 并取强度
        [wS0, ~, ~, ~] = cumulativeQUV(IMG1_wholeStr, IMG2_wholeStr);
        strLin = mean(wS0, 3);
        SS = 20 * log10(strLin + eps);

        global_min = min(global_min, min(SS(:)));
        global_max = max(global_max, max(SS(:)));
    end

    if isinf(global_min) || isinf(global_max)
        warning('无法获得有效的强度范围，跳过文件: %s', infile);
        fclose(fid);
        continue;
    end

    strLrg = global_min + 5;
    strUrg = global_max - 5;
    if strUrg <= strLrg
        % 避免除0
        strLrg = global_min;
        strUrg = global_max;
    end

    fprintf('归一化范围: [%.3f, %.3f]\n', strLrg, strUrg);

    % 第二遍：使用并行（parfor）并行生成每帧 PNG
    fprintf('\n第二遍扫描 (保存 PNG)...\n');

    frame_idx = 1:nY;

    if use_parallel
        parfor_progress = 0; % 仅用于简单计数（不能直接打印过多）
        parfor fi_idx = 1:length(frame_idx)
            iy = frame_idx(fi_idx);
            try
                fidp = fopen(infile, 'r');
                if fidp==-1, error('无法打开文件在 worker'); end
                fseek(fidp, bob + (SPL * nX * 2 + 2) * 2 * (nR * (iy - 1)), 'bof');
                Bs1 = zeros(SPL, nX, nR, 'single');
                Bs2 = zeros(SPL, nX, nR, 'single');
                for Ic = 1:nR
                    fseek(fidp, 4, 'cof');
                    B = single(reshape(fread(fidp, SPL * nX * 2, 'int16'), [SPL, nX * 2]));
                    Bs1(:, :, Ic) = (B(:, 1:nX) - single(Ref_ch1));
                    Bs2(:, :, Ic) = (B(:, nX + 1:end) - single(Ref_ch2));
                end
                fclose(fidp);

                if use_gpu
                    % 可选：使用 gpu 加速 FFT（注意 cumulativeQUV 可能不支持 gpuArray）
                    try
                        gBs1 = gpuArray(Bs1);
                        gBs2 = gpuArray(Bs2);
                        gBimg1 = fft(gBs1 .* reshape(gpuArray(winG_whole), [SPL,1,1]), SPL, 1);
                        gBimg2 = fft(gBs2 .* reshape(gpuArray(winG_whole), [SPL,1,1]), SPL, 1);
                        IMG1_wholeStr = gather(gBimg1(czrg, :, :));
                        IMG2_wholeStr = gather(gBimg2(czrg, :, :));
                    catch
                        % 回退到 CPU 路径
                        Bimg1_wholeStr = fft(Bs1 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
                        Bimg2_wholeStr = fft(Bs2 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
                        IMG1_wholeStr = Bimg1_wholeStr(czrg, :, :);
                        IMG2_wholeStr = Bimg2_wholeStr(czrg, :, :);
                    end
                else
                    Bimg1_wholeStr = fft(Bs1 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
                    Bimg2_wholeStr = fft(Bs2 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
                    IMG1_wholeStr = Bimg1_wholeStr(czrg, :, :);
                    IMG2_wholeStr = Bimg2_wholeStr(czrg, :, :);
                end

                [wS0, ~, ~, ~] = cumulativeQUV(IMG1_wholeStr, IMG2_wholeStr);
                strLin = mean(wS0, 3);
                SS = 20 * log10(strLin + eps);

                % 归一化并裁剪
                Struc_loc = (SS - strLrg) ./ (strUrg - strLrg);
                Struc_loc = min(max(Struc_loc, 0), 1);

                % 保存 PNG
                out_fname = fullfile(outdir, sprintf('%s_frame_%04d.png', name, iy));
                try
                    imwrite(uint8(255 * Struc_loc), out_fname);
                catch ME
                    warning('保存 PNG 失败: %s, 错误: %s', out_fname, ME.message);
                end
            catch ME
                warning('帧 %d 处理失败: %s', iy, ME.message);
            end
        end

    else
        for iy = 1:nY
            if mod(iy, max(1, floor(nY/10))) == 0 || iy == nY
                fprintf('  处理帧 %d/%d\n', iy, nY);
            end
            fseek(fid, bob + (SPL * nX * 2 + 2) * 2 * (nR * (iy - 1)), 'bof');
            Bs1 = zeros(SPL, nX, nR, 'single');
            Bs2 = zeros(SPL, nX, nR, 'single');
            for Ic = 1:nR
                fseek(fid, 4, 'cof');
                B = single(reshape(fread(fid, SPL * nX * 2, 'int16'), [SPL, nX * 2]));
                Bs1(:, :, Ic) = (B(:, 1:nX) - single(Ref_ch1));
                Bs2(:, :, Ic) = (B(:, nX + 1:end) - single(Ref_ch2));
            end

            Bimg1_wholeStr = fft(Bs1 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
            Bimg2_wholeStr = fft(Bs2 .* reshape(winG_whole, [SPL,1,1]), SPL, 1);
            IMG1_wholeStr = Bimg1_wholeStr(czrg, :, :);
            IMG2_wholeStr = Bimg2_wholeStr(czrg, :, :);

            [wS0, ~, ~, ~] = cumulativeQUV(IMG1_wholeStr, IMG2_wholeStr);
            strLin = mean(wS0, 3);
            SS = 20 * log10(strLin + eps);

            % 归一化并裁剪
            Struc = (SS - strLrg) ./ (strUrg - strLrg);
            Struc = min(max(Struc, 0), 1);

            % 保存 PNG (Z x X) -> 转换为 uint8
            out_fname = fullfile(outdir, sprintf('%s_frame_%04d.png', name, iy));
            try
                imwrite(uint8(255 * Struc), out_fname);
            catch ME
                warning('保存 PNG 失败: %s, 错误: %s', out_fname, ME.message);
            end
        end
    end

    fclose(fid);
    fprintf('完成 %s，PNG 已保存到 %s\n', name, outdir);
end

fprintf('\n所有文件处理完成！\n');
