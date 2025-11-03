function topLines = load_topLines_mat(matfile, nX, nY, nZ)
% load_topLines_mat 从 .mat 文件加载 topLines（或同名变量），并调整为 [nX x nY] 尺寸，且深度在 [1,nZ]
%   topLines = load_topLines_mat(matfile, nX, nY, nZ)
%   支持多种变量名和尺寸适配策略（transpose、reshape、imresize、interp2）

    if ~exist(matfile, 'file')
        error('指定的 mat 文件不存在: %s', matfile);
    end
    S = load(matfile);

    % 可能的变量名候选
    cand = {'topLines','TopLines','toplines','surface','surface_line','surfaceLine'};
    varName = '';
    fn = fieldnames(S);
    for i=1:numel(cand)
        if ismember(cand{i}, fn)
            varName = cand{i}; break;
        end
    end
    if isempty(varName)
        % 如果只有一个变量，尝试使用它
        if numel(fn) == 1
            varName = fn{1};
        else
            error('未在 %s 中找到 topLines 变量，请使用变量名 topLines 或指定 matSegFile', matfile);
        end
    end

    val = S.(varName);
    % 如果是结构体，尝试找到内部字段
    if isstruct(val)
        f = fieldnames(val);
        if ismember('topLines', f)
            val = val.topLines;
        elseif ismember('TopLines', f)
            val = val.TopLines;
        else
            error('在 %s 中找到结构变量 %s，但未包含 topLines 字段', matfile, varName);
        end
    end

    % 转为双精度矩阵
    val = double(val);

    % 处理不同尺寸情况，目标是得到 [nX x nY]
    sz = size(val);
    if isequal(sz, [nX, nY])
        topLines = val;
    elseif isequal(sz, [nY, nX])
        topLines = val';
    elseif isvector(val)
        v = val(:)';
        if numel(v) == nX
            % 单列的 topLine 向量，复制到所有 B-scan
            topLines = repmat(v', 1, nY);
        elseif numel(v) == nY
            % 单行，为每个B-scan 指定相同深度(不常见)
            topLines = repmat(v, nX, 1);
        elseif numel(v) == nX*nY
            topLines = reshape(v, [nX, nY]);
        else
            error('加载的向量长度 (%d) 无法匹配目标尺寸 [%d x %d]', numel(v), nX, nY);
        end
    else
        % 如果列数与 nY 不同，则在列方向上插值/重采样
        if sz(1) == nX && sz(2) ~= nY
            % 沿列插值
            topLines = zeros(nX, nY);
            oldX = 1:sz(2);
            newX = linspace(1, sz(2), nY);
            for r = 1:nX
                row = val(r, :);
                topLines(r, :) = interp1(oldX, row, newX, 'linear', 'extrap');
            end
        elseif sz(2) == nY && sz(1) ~= nX
            % 沿行插值
            oldY = 1:sz(1);
            newY = linspace(1, sz(1), nX);
            tmp = zeros(nX, nY);
            for c = 1:nY
                col = val(:, c);
                tmp(:, c) = interp1(oldY, col, newY, 'linear', 'extrap');
            end
            topLines = tmp;
        else
            % 通用重采样：使用 imresize 如果可用，否则双线性插值
            try
                topLines = imresize(val, [nX, nY], 'bilinear');
            catch
                % 手动双线性插值
                [r_old, c_old] = size(val);
                [C, R] = meshgrid(1:c_old, 1:r_old);
                [Cq, Rq] = meshgrid(linspace(1,c_old,nY), linspace(1,r_old,nX));
                topLines = interp2(C, R, val, Cq, Rq, 'linear');
            end
        end
    end

    % 限制并取整到有效范围
    topLines = round(topLines);
    topLines(isnan(topLines)) = 1;
    topLines = max(topLines, 1);
    topLines = min(topLines, nZ);
end
