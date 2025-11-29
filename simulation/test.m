close all; clear; clc;

%% ==============================
%  通用参数
%% ==============================
saveDir = fullfile(pwd, 'output');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

rng(1); % 固定随机种子

%% =============================================================
%  函数：绘制单位球 + 点 + 切圆
%% =============================================================
function drawOneCase(P, figName, saveDir)
    figure('Color','w','Position',[200 200 600 600]); hold on; axis equal;

    % 绘制单位球
    [xs,ys,zs] = sphere(80);
    surf(xs,ys,zs,'FaceAlpha',0.05,'EdgeColor','none'); 
    colormap(gray);

    % 绘制三个点
    scatter3(P(:,1),P(:,2),P(:,3),120,'r','filled','MarkerEdgeColor','k');

    % 拟合平面
    [n,c] = planeFit(P);

    % 画切圆
    drawCutCircle(n,c,200);

    xlabel('Q'); ylabel('U'); zlabel('V');
    title(figName);
    view(3); axis vis3d;

    % 保存图像
    outPNG = fullfile(saveDir, [figName '.png']);
    saveas(gcf, outPNG);
    fprintf("Saved: %s\n", outPNG);
end


%% =============================================================
%  平面拟合函数：n·x = c
%% =============================================================
close all; clear; clc;

%% ==============================
%  科研风格庞加莱球可视化（主脚本）
%% =============================
saveDir = fullfile(pwd, 'output');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

rng(1); % 固定随机种子

%% =============================================================
% Case 1：强双折射（三个点相距较远）- 不同方向
%% =============================================================
theta = [0.6, 0.9, 1.3];  % 更换角度
phi   = [1.5, 3.0, 4.8];  % 更换方位
P1 = [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)]';

drawPoincareSphere(P1, 'Strong_Birefringence', saveDir);

%% =============================================================
% Case 2：微弱双折射（三个点较接近但可见）- 不同方向
%% =============================================================
center = [0.8 0.3 0.5]; center = center/norm(center);  % 侧面位置
eps_small = 0.06;  % 较小的扰动
P2 = repmat(center, 3, 1) + eps_small*randn(3, 3);
P2 = P2 ./ vecnorm(P2, 2, 2);

drawPoincareSphere(P2, 'Weak_Birefringence', saveDir);

%% =============================================================
% Case 3：对比演示 - 正常拟合 vs 强制大半径拟合
%% =============================================================
fprintf('\n========== 对比演示：强制大半径的问题 =========\n');

% 使用弱双折射的数据
P_demo = P2;

% 正常拟合
[n_normal, c_normal, r_normal, error_normal] = planeFit(P_demo);
fprintf('【正常拟合】\n');
fprintf('  切圆半径 r = %.4f\n', r_normal);
fprintf('  拟合误差 = %.6f\n', error_normal);

% 强制半径 >= 0.8（接近大圆）
r_min = 0.8;
[n_forced, c_forced, r_forced, error_forced] = planeFitConstrained(P_demo, r_min);
fprintf('\n【强制半径 >= %.2f】\n', r_min);
fprintf('  切圆半径 r = %.4f\n', r_forced);
fprintf('  拟合误差 = %.6f\n', error_forced);
fprintf('  误差增大了 %.1f 倍！\n', error_forced / max(error_normal, eps));

% 绘制对比图
drawComparisonFigure(P_demo, n_normal, c_normal, r_normal, error_normal, ...
                      n_forced, c_forced, r_forced, error_forced, r_min, saveDir);

%% =============================================================
% 以下为本脚本使用的本地函数（放在脚本末尾，兼容 MATLAB 脚本风格）
%% =============================================================
function drawPoincareSphere(P, figName, saveDir)
    figure('Color','w','Position',[200 200 800 800]);
    hold on; axis equal;

    % 1. 透明球壳
    [sxs, sys, szs] = sphere(80);
    surf(sxs, sys, szs, 'FaceColor', [0.96 0.97 1], 'FaceAlpha', 0.06, 'EdgeColor', 'none');

    % 2. 十字坐标系虚线
    plot3([-1.1 1.1], [0 0], [0 0], 'k--', 'LineWidth', 1.2);
    plot3([0 0], [-1.1 1.1], [0 0], 'k--', 'LineWidth', 1.2);
    plot3([0 0], [0 0], [-1.1 1.1], 'k--', 'LineWidth', 1.2);

    % 3. 两条主要大圈（Q-U, Q-V）
    t = linspace(0, 2*pi, 300);
    plot3(cos(t), sin(t), zeros(size(t)), 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
    plot3(cos(t), zeros(size(t)), sin(t), 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);

    % 标签
    text(1.08, 0, -1.12, 'Q', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, 1.08, -1.12, 'U', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, 0, 1.18, 'V', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');

    % 拟合平面并计算法向量
    [n, c, r, error] = planeFit(P);
    fprintf('  切圆半径: r = %.4f, 拟合误差: %.6f\n', r, error);

    % 绘制切圆和平面
    drawCutCircleAndPlane(n, c, 200);

    % 绘制点
    scatter3(P(:,1), P(:,2), P(:,3), 200, [0.9 0.1 0.1], 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);

    % 视角和光照
    view([-40, 25]);
    camlight('left'); camlight('right');
    lighting phong;
    material([0.4 0.6 0.5 10 0.8]);

    axis([-1 1 -1 1 -1 1]);
    set(gca, 'Projection', 'perspective');
    grid on; box on;
    set(gca, 'XTick', [-1 -0.5 0 0.5 1], 'YTick', [-1 -0.5 0 0.5 1], 'ZTick', [-1 -0.5 0 0.5 1]);
    set(gca, 'FontSize', 10);

    title(figName, 'FontSize', 16, 'FontWeight', 'bold');

    outPNG = fullfile(saveDir, [figName '.png']);
    saveas(gcf, outPNG);
    fprintf('✅ Saved: %s\n', outPNG);
end

function [n, c, r, error] = planeFit(P)
    Pm = mean(P, 1);
    Pc = P - Pm;
    [~, ~, V] = svd(Pc, 0);
    n = V(:, end);
    n = n / norm(n);
    c = dot(n, Pm);
    % clamp c to [-1,1] for numerical stability
    c = max(min(c, 1), -1);
    r = sqrt(max(0, 1 - c^2));
    error = sum((P * n - c).^2);
end

function [n, c, r, error] = planeFitConstrained(P, r_min)
    % 约束：|c| <= sqrt(1 - r_min^2)
    c_max = sqrt(max(0, 1 - r_min^2));

    % 初始猜测
    [n0, c0, ~, err0] = planeFit(P);

    % 如果初始解已经满足约束，直接返回
    if sqrt(max(0, 1 - c0^2)) >= r_min
        n = n0; c = c0; r = sqrt(max(0, 1 - c^2)); error = err0; return;
    end

    % 优化变量 x = [n(1:3); c]
    objective = @(x) sum((P * x(1:3) - x(4)).^2);

    % 非线性约束：norm(n)==1, |c|<=c_max -> represented as inequalities
    function [cineq, ceq]
        % empty - we use bounds for c and a normalization equality
    end

    % Use simple parameterization: enforce unit length by normalizing inside objective
    wrappedObj = @(x) objective([x(1:3)/norm(x(1:3)); x(4)]);

    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxIterations', 500);
    x0 = [n0; c0];
    lb = [-inf; -inf; -inf; -c_max];
    ub = [inf; inf; inf; c_max];

    x_opt = fmincon(wrappedObj, x0, [], [], [], [], lb, ub, [], options);

    n = x_opt(1:3); n = n / norm(n);
    c = x_opt(4);
    c = max(min(c, 1), -1);
    r = sqrt(max(0, 1 - c^2));
    error = sum((P * n - c).^2);
end

function drawCutCircleAndPlane(n, c, N)
    if abs(c) > 1
        warning('Plane does not intersect unit sphere.');
        return;
    end
    r = sqrt(max(0, 1 - c^2));

    if abs(n(1)) < 0.9
        a = [1; 0; 0];
    else
        a = [0; 1; 0];
    end
    u = a - n*(n'*a); u = u / norm(u);
    v = cross(n, u);

    t = linspace(0, 2*pi, N);
    circle = c*n + r*(u*cos(t) + v*sin(t));

    % 填充圆盘
    theta = linspace(0, 2*pi, 60);
    rho = linspace(0, r, 30);
    [T, R] = meshgrid(theta, rho);
    X = c*n(1) + R.*(u(1)*cos(T) + v(1)*sin(T));
    Y = c*n(2) + R.*(u(2)*cos(T) + v(2)*sin(T));
    Z = c*n(3) + R.*(u(3)*cos(T) + v(3)*sin(T));
    surf(X, Y, Z, 'FaceColor', [0.98 0.93 0.6], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    % 切圆边界
    plot3(circle(1,:), circle(2,:), circle(3,:), 'Color', [0.85 0.45 0.1], 'LineWidth', 1.5);
    idx = 1:2:length(t);
    plot3(circle(1,idx), circle(2,idx), circle(3,idx), ':', 'Color', [0.85 0.45 0.1], 'LineWidth', 1);

    % 法向量箭头
    center = c * n;
    arrowLen = 0.45;
    quiver3(center(1), center(2), center(3), n(1)*arrowLen, n(2)*arrowLen, n(3)*arrowLen, ...
        'Color', [0.8 0.2 0.2], 'LineWidth', 1.2, 'MaxHeadSize', 0.7, 'AutoScale', 'off');
end

function drawComparisonFigure(P, n1, c1, r1, err1, n2, c2, r2, err2, r_min, saveDir)
    figure('Color','w','Position',[100 100 1600 700]);

    % 左图：正常拟合
    subplot(1,2,1); hold on; axis equal;
    [sxsL, sysL, szsL] = sphere(80);
    surf(sxsL, sysL, szsL, 'FaceColor', [0.96 0.97 1], 'FaceAlpha', 0.06, 'EdgeColor', 'none');
    plot3([-1.1 1.1], [0 0], [0 0], 'k--', 'LineWidth', 1.2);
    plot3([0 0], [-1.1 1.1], [0 0], 'k--', 'LineWidth', 1.2);
    plot3([0 0], [0 0], [-1.1 1.1], 'k--', 'LineWidth', 1.2);
    t = linspace(0, 2*pi, 300);
    plot3(cos(t), sin(t), zeros(size(t)), 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
    plot3(cos(t), zeros(size(t)), sin(t), 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
    drawCutCircleAndPlane(n1, c1, 150);
    text(1.08, 0, -1.12, 'Q', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, 1.08, -1.12, 'U', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, 0, 1.18, 'V', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    center1 = c1 * n1;
    scatter3(center1(1), center1(2), center1(3), 150, [0 0.6 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    scatter3(P(:,1), P(:,2), P(:,3), 250, [0.9 0.1 0.1], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    for i = 1:size(P, 1)
        proj = P(i,:)' - n1 * (dot(P(i,:), n1) - c1);
        plot3([P(i,1) proj(1)], [P(i,2) proj(2)], [P(i,3) proj(3)], 'g--', 'LineWidth', 1.5);
    end
    view([30, 20]); camlight('left'); camlight('right'); lighting phong;
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]); grid on; box on; set(gca, 'FontSize', 11);
    title(sprintf('正常拟合\n半径 r=%.3f, 误差=%.6f', r1, err1), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0.6 0]);
    xlabel('Q', 'FontSize', 12, 'FontWeight', 'bold'); ylabel('U', 'FontSize', 12, 'FontWeight', 'bold'); zlabel('V', 'FontSize', 12, 'FontWeight', 'bold');

    % 右图：强制大半径拟合
    subplot(1,2,2); hold on; axis equal;
    [sxsR, sysR, szsR] = sphere(80);
    surf(sxsR, sysR, szsR, 'FaceColor', [0.96 0.97 1], 'FaceAlpha', 0.06, 'EdgeColor', 'none');
    plot3([-1.1 1.1], [0 0], [0 0], 'k--', 'LineWidth', 1.2);
    plot3([0 0], [-1.1 1.1], [0 0], 'k--', 'LineWidth', 1.2);
    plot3([0 0], [0 0], [-1.1 1.1], 'k--', 'LineWidth', 1.2);
    t2 = linspace(0, 2*pi, 300);
    plot3(cos(t2), sin(t2), zeros(size(t2)), 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
    plot3(cos(t2), zeros(size(t2)), sin(t2), 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
    text(1.08, 0, -1.12, 'Q', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, 1.08, -1.12, 'U', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, 0, 1.18, 'V', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
    drawCutCircleAndPlane(n2, c2, 150);
    center2 = c2 * n2;
    scatter3(center2(1), center2(2), center2(3), 150, [0.8 0 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    center1 = c1 * n1;
    plot3([center1(1) center2(1)], [center1(2) center2(2)], [center1(3) center2(3)], 'Color', [0.5 0 0.5], 'LineWidth', 2.5, 'LineStyle', '--');
    scatter3(P(:,1), P(:,2), P(:,3), 250, [0.9 0.1 0.1], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    for i = 1:size(P, 1)
        proj = P(i,:)' - n2 * (dot(P(i,:), n2) - c2);
        plot3([P(i,1) proj(1)], [P(i,2) proj(2)], [P(i,3) proj(3)], 'r--', 'LineWidth', 2);
    end
    view([30, 20]); camlight('left'); camlight('right'); lighting phong;
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]); grid on; box on; set(gca, 'FontSize', 11);

    error_ratio = err2 / max(err1, 1e-10);
    title(sprintf('强制半径≥%.2f\n半径 r=%.3f, 误差=%.6f\n(误差增大 %.0f 倍)', r_min, r2, err2, error_ratio), ...
          'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
    xlabel('Q', 'FontSize', 12, 'FontWeight', 'bold'); ylabel('U', 'FontSize', 12, 'FontWeight', 'bold'); zlabel('V', 'FontSize', 12, 'FontWeight', 'bold');

    outPNG = fullfile(saveDir, 'Comparison_Normal_vs_Forced.png');
    saveas(gcf, outPNG);
    fprintf('\n✅ 对比图已保存: %s\n', outPNG);
end
