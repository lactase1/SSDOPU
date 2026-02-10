% extract_dcm_frame.m
% 提取指定帧并把符合名字规则的 DCM 转为 TIFF 和 PNG
% 直接在脚本中指定路径和目标帧，运行时无需传参

% --- 配置（按需修改） ---
dcm_folder = 'D:\1-Liu Jian\yongxin.wang\Output\Enface\ddg_3layer_3_3\2024.09.04_18.42.53_1.SLnoBG_disk\dcm';
target_frame = 26; % 1-based
output_dir = 'D:\1-Liu Jian\yongxin.wang\Output\Enface\ddg_3layer_3_3\2024.09.04_18.42.53_1.SLnoBG_disk\slect_frame';
% -------------------------

if ~exist(dcm_folder, 'dir')
    error('DCM 文件夹不存在: %s', dcm_folder);
end

tiff_dir = fullfile(output_dir, 'tiff');
png_dir  = fullfile(output_dir, 'png');
if ~exist(tiff_dir, 'dir'), mkdir(tiff_dir); end
if ~exist(png_dir, 'dir'), mkdir(png_dir); end

files = dir(fullfile(dcm_folder, '*.dcm'));
fprintf('查找 DCM：%d 个文件，开始按规则筛选...\n', numel(files));

% 过滤：包含 '2-2'、'2-5'（但排除含 'flat' 的展平文件），或包含 '_3-'（匹配 3-x 形式，如 _3-3_）
selected = {};
for k = 1:numel(files)
    nm = files(k).name;
    nm_low = lower(nm);
    is1_1 = contains(nm, '1-1') && ~contains(nm_low, 'flat');
    is2_2 = contains(nm, '2-2') && ~contains(nm_low, 'flat');
    is2_5 = contains(nm, '2-5') && ~contains(nm_low, 'flat');
    is3 = contains(nm, '_3-') || contains(nm, '-3-') || contains(nm, '-3_');
    is4 = contains(nm, '_4-') || contains(nm, '-4-') || contains(nm, '-4_');
    if is1_1 || is2_2 || is2_5 || is3 || is4
        selected{end+1} = fullfile(dcm_folder, nm); %#ok<AGROW>
    end
end

fprintf('筛选后 %d 个文件，将提取第 %d 帧并保存到 %s (TIFF) 和 %s (PNG)\n', numel(selected), target_frame, tiff_dir, png_dir);

success = 0;
fail_list = {};
for i = 1:numel(selected)
    fpath = selected{i};
    [~, base, ~] = fileparts(fpath);
    try
        % 优先使用 dicomread 直接按帧读取（若支持）
        info = dicominfo(fpath);
        % 计算帧索引
        if isfield(info, 'NumberOfFrames') && info.NumberOfFrames > 1
            frame_idx = min(target_frame, info.NumberOfFrames);
            try
                img = dicomread(info, 'frames', frame_idx); % 按帧读取（MATLAB 版本支持时）
            catch
                % 如果按帧读取失败，尝试读取全部再切片
                fullimg = dicomread(info);
                % fullimg 可能是 3D 或 4D
                if ndims(fullimg) == 3
                    img = fullimg(:, :, frame_idx);
                elseif ndims(fullimg) == 4
                    img = squeeze(fullimg(:, :, :, frame_idx));
                else
                    img = squeeze(fullimg(:,:,frame_idx));
                end
            end
        else
            % 单帧 DCM
            img = dicomread(info);
        end
    catch ME
        % 回退：尝试手工解析 PixelData（针对不被 dicomread 支持的 bit-packed 数据，例如 1-bit）
        try
            info = dicominfo(fpath);
            bits = double(getfield(info, 'BitsAllocated'));
            rows = double(getfield(info, 'Rows'));
            cols = double(getfield(info, 'Columns'));
            if isfield(info, 'NumberOfFrames')
                frames = double(info.NumberOfFrames);
            else
                frames = 1;
            end

            raw = fopen(fpath, 'r'); fclose(raw); % 确保文件可访问
            % 注意：MATLAB 的 dicominfo 没有直接提供 PixelData 字节流，这里仍尝试用 dicomread 作为最后手段
            img_all = dicomread(info); % 再一次尝试 (有时 dicomread 在 info 形式会成功)
            if ndims(img_all) == 3
                frame_idx = min(target_frame, size(img_all,3));
                img = img_all(:,:,frame_idx);
            elseif ndims(img_all) == 4
                frame_idx = min(target_frame, size(img_all,4));
                img = squeeze(img_all(:,:,:,frame_idx));
            else
                img = img_all;
            end
        catch ME2
            fprintf('处理文件 %s 出错 (读取像素失败): %s\n', files(i).name, ME2.message);
            fail_list{end+1} = files(i).name; %#ok<AGROW>
            continue;
        end
    end

    % 把 img 转为可写入的 uint8/uint8 RGB
    try
        % 如果是逻辑型或二值图
        if islogical(img)
            img_u8 = uint8(img) * 255;
        elseif ndims(img) == 3 && size(img,3) == 3 && isa(img, 'uint8')
            % RGB uint8，直接使用
            img_u8 = img;
        else
            % 灰度或其他位深，缩放到 0-255
            img_double = double(img);
            mn = min(img_double(:));
            mx = max(img_double(:));
            if mx == mn
                % 常量值：扩展为小图像（避免图像太小）
                v = uint8(0);
                if ~isnan(mn)
                    v = uint8( round( (mn - mn) / max(1, mx-mn) * 255 ) );
                end
                img_u8 = repmat(v, [64, 64]);
            else
                img_u8 = uint8( round( (img_double - mn) / (mx - mn) * 255 ) );
            end
        end
    catch ME
        fprintf('处理像素并归一化时出错: %s\n', ME.message);
        fail_list{end+1} = files(i).name; %#ok<AGROW>
        continue;
    end

    % 如果图像仍然是 1x1（单像素），扩展为 64x64
    if ismatrix(img_u8) && all(size(img_u8) == [1 1])
        img_u8 = repmat(img_u8, 64, 64);
    end

    % 输出路径
    tiff_name = sprintf('%s_frame%d.tiff', base, target_frame);
    png_name  = sprintf('%s_frame%d.png', base, target_frame);

    try
        imwrite(img_u8, fullfile(tiff_dir, tiff_name));
        imwrite(img_u8, fullfile(png_dir, png_name));
        success = success + 1;
    catch ME
        fprintf('保存 %s/%s 失败: %s\n', tiff_name, png_name, ME.message);
        fail_list{end+1} = files(i).name; %#ok<AGROW>
        continue;
    end
end

% 精简输出（两行）
fprintf('\n转换完成：成功 %d / %d 个文件。\n', success, numel(selected));
if ~isempty(fail_list)
    fprintf('失败文件（%d 个）：%s\n', numel(fail_list), strjoin(fail_list, ', '));
end
fprintf('输出目录：%s (tiff) 和 %s (png)\n', tiff_dir, png_dir);

% End of script
