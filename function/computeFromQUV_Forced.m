function [opticAxis, phaseRetardation] = computeFromQUV_Forced(QUV, params)
% 使用强制大半径约束法计算光轴和相位延迟
% QUV: [depth x width x 3]
% 输出: opticAxis [depth_out x width x 3], phaseRetardation [depth_out-1 x width]

[depth_in, width, ~] = size(QUV);
r_min = params.r_min;
depth_out = params.output_depth_forced;
window_size = params.window_forced_max;
half_win = floor(window_size / 2);

opticAxis = zeros(depth_out, width, 3);
phaseRetardation = zeros(depth_out - 1, width);

centers_f = linspace(1, depth_in, depth_out);
reference_B = [];

for w = 1:width
    points = squeeze(QUV(:, w, :));

    for d = 1:depth_out
        center_idx = round(centers_f(d));
        start_idx = max(1, center_idx - half_win);
        end_idx = start_idx + window_size - 1;
        if end_idx > depth_in
            end_idx = depth_in;
            start_idx = max(1, end_idx - window_size + 1);
        end

        window_points = points(start_idx:end_idx, :);
        [B, ~, ~, ~] = planeFitConstrained_local(window_points, r_min);

        if isempty(reference_B)
            reference_B = B;
        elseif dot(B, reference_B) < 0
            B = -B;
        end

        opticAxis(d, w, :) = B;
    end

    for d = 1:depth_out - 1
        center_idx = round(centers_f(d));
        start_idx = max(1, center_idx - half_win);
        end_idx = start_idx + window_size - 1;
        if end_idx > depth_in
            end_idx = depth_in;
            start_idx = max(1, end_idx - window_size + 1);
        end

        if start_idx < depth_in
            P1 = points(start_idx, :)';
            P2 = points(min(start_idx + 1, depth_in), :)';
            dot_product = dot(P1, P2);
            dot_product = max(-1, min(1, dot_product));
            phaseRetardation(d, w) = 0.5 * acosd(dot_product);
        else
            phaseRetardation(d, w) = 0;
        end
    end
end
end

function [n, c, r, error] = planeFitConstrained_local(P, r_min)
% 拟合切平面并强制半径不小于 r_min（球半径为1）
    if isempty(P) || size(P, 1) < 3
        n = [0; 0; 1]; c = 0; r = 1; error = 0; return;
    end

    [n0, c0, ~, err0] = planeFit_local(P);
    r0 = sqrt(max(0, 1 - c0^2));
    c_max = sqrt(max(0, 1 - r_min^2));

    if r0 >= r_min
        n = n0; c = c0; r = r0; error = err0; return;
    end

    c = sign(c0) * c_max;
    n = n0;
    r = sqrt(max(0, 1 - c^2));
    error = sum((P * n - c).^2);
end

function [n, c, r, error] = planeFit_local(P)
% 基础切平面拟合（SVD）
    Pm = mean(P, 1);
    Pc = P - Pm;
    [~, ~, V] = svd(Pc, 0);
    n = V(:, end);
    n = n / norm(n);
    c = dot(n, Pm);
    c = max(min(c, 1), -1);
    r = sqrt(max(0, 1 - c^2));
    error = sum((P * n - c).^2);
end
