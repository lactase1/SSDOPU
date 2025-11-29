%% fig1.m: 基于Polarization State Tracing (PST) 方法模拟庞加莱球轨迹

% 清理环境
close all; clear; clc;

%% 定义Mueller矩阵函数
function M = mueller_matrix(axisVec, delta)
    % Mueller 矩阵（理想无耗散延迟器）
    % Retarder 在 Q-U-V 子空间等价为绕某轴的 3D 旋转，I 分量不变。
    % axisVec: 3x1 单位向量指定旋转轴在 QUV 子空间的位置
    a = axisVec(:);
    a = a / (norm(a) + eps);
    % Rodrigues 旋转矩阵（3x3）
    K = [   0, -a(3),  a(2);
         a(3),    0, -a(1);
        -a(2), a(1),    0];
    R = cos(delta) * eye(3) + (1 - cos(delta)) * (a * a') + sin(delta) * K;
    M = [1, 0, 0, 0;
         0, R(1,1), R(1,2), R(1,3);
         0, R(2,1), R(2,2), R(2,3);
         0, R(3,1), R(3,2), R(3,3)];
end

%% 绘制庞加莱球
figure;
[X, Y, Z] = sphere(50); % 创建单位球面数据
surf(X, Y, Z, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1); % 绘制球面
axis equal;
xlabel('Q'); ylabel('U'); zlabel('V');
grid on; hold on;

%% 模拟层结构
numLayers = 6; % 层数
% 光轴：使用 3D 向量（Q,U,V）——选择不共线的两个方向以演示 PST
opticAxes = [
    0.7071, 0.0000, 0.7071;  % axis1 (Q-V 混合)
    0.7071, 0.0000, 0.7071;
    0.7071, 0.0000, 0.7071;
    0.0000, 0.7071, 0.7071;  % axis2 (U-V 混合)
    0.0000, 0.7071, 0.7071;
    0.0000, 0.7071, 0.7071;
];
phaseRetardations = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]; % 每层的相位延迟

% 初始输入斯托克斯向量（必须为完全偏振态，避免落到原点导致除零）
% 选择一个不与光轴共线的完全偏振态：Q=0.5, U=0.5, V=0.70710678 (Q^2+U^2+V^2 = 1)
S_in = [1; 0.5; 0.5; 0.70710678];

% 计算输出斯托克斯参数（累积变换）
stokesOutputs = zeros(numLayers+1, 4);
stokesOutputs(1, :) = S_in'; % P0: 表面

M_cum = eye(4);
for i = 1:numLayers
    axisVec = opticAxes(i, :)';
    delta = phaseRetardations(i);
    M = mueller_matrix(axisVec, delta);
    M_cum = M * M_cum;
    S_out = M_cum * S_in;
    stokesOutputs(i+1, :) = S_out';
end

% 计算庞加莱球坐标 (Q/I, U/I, V/I)
poincarePoints = stokesOutputs ./ stokesOutputs(:, 1);

% 绘制斯托克斯参数轨迹（更清晰的可视化）
trajH = plot3(poincarePoints(:, 2), poincarePoints(:, 3), poincarePoints(:, 4), '-','Color',[0.85 0.2 0.2],'LineWidth',1.5);
hold on;
% 标注每个采样点（P0..Pn）并用不同颜色和大小突出
colors = lines(size(poincarePoints,1));
for k = 1:size(poincarePoints,1)
    scatter3(poincarePoints(k,2), poincarePoints(k,3), poincarePoints(k,4), 60, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor','k');
    text(poincarePoints(k,2)+0.02, poincarePoints(k,3)+0.02, poincarePoints(k,4)+0.02, sprintf('P_%d', k-1), 'FontSize', 9);
end

% 把球面透明度调小并确保轨迹点在前景
uistack(trajH,'top');

%% 计算局部光轴和相位延迟使用PST方法
% 假设每3层一组相同光轴，前3层theta=0，后3层pi/2

% 拟合平面法向量函数
function normal = fit_plane_normal(p1, p2, p3)
    v1 = p2 - p1;
    v2 = p3 - p1;
    n = cross(v1, v2);
    nnorm = norm(n);
    if nnorm < 1e-8
        % 点近共线或重复：退化到使用 SVD / PCA 风格的近似。
        pts = [p1(:)'; p2(:)'; p3(:)'];
        % 中心化
        pts0 = pts - mean(pts,1);
        % SVD 找到最小特征向量作为法向量
        [~, ~, V] = svd(pts0, 0);
        n = V(:,end)';
        nnorm = norm(n);
        if nnorm < 1e-8
            % 最后兜底，返回 z 轴
            normal = [0 0 1];
            return;
        end
    end
    normal = n / nnorm;
end

% 旋转向量函数 (Rodrigues)
function v_rot = rotate_vector(v, axis, angle)
    axis = axis / norm(axis);
    v_rot = v * cos(angle) + cross(axis, v) * sin(angle) + axis * dot(axis, v) * (1 - cos(angle));
end

% 计算局部相位延迟函数 (Eq 3.5)
function delta = calculate_local_retardation(P_prev, P_curr, A_prime)
    % 计算由 Eq.(3.5) 得到的半角值：delta = 0.5 * acos( (cross(P_prev,A) · cross(P_curr,A)) / (|cross(P_prev,A)| |cross(P_curr,A)|) )
    cross_prev = cross(P_prev, A_prime);
    cross_curr = cross(P_curr, A_prime);
    n1 = norm(cross_prev);
    n2 = norm(cross_curr);
    epsval = 1e-10;
    if n1 < epsval || n2 < epsval
        % 如果任意 cross 近零，说明 P 与 A_prime 共线或近共线，半角退化为0
        delta = 0;
        return;
    end
    cos_angle = dot(cross_prev, cross_curr) / (n1 * n2);
    % 数值稳定：夹紧到 [-1,1]
    cos_angle = max(-1, min(1, cos_angle));
    delta = 0.5 * acos(cos_angle);
end

    % 绘制一个半透明平面的辅助函数
    function h = plot_plane(center, normal, color)
        % 将 normal 单位化
        normal = normal / (norm(normal) + eps);
        % 构造局部平面坐标系（找两个正交切向量）
        if abs(normal(3)) < 0.9
            tangent1 = cross(normal, [0 0 1]);
        else
            tangent1 = cross(normal, [0 1 0]);
        end
        tangent1 = tangent1 / norm(tangent1 + eps);
        tangent2 = cross(normal, tangent1);
        tangent2 = tangent2 / (norm(tangent2) + eps);
        % 生成平面网格
        [u,v] = meshgrid(linspace(-0.35,0.35,20));
        pts = center(:)' + u(:)*tangent1 + v(:)*tangent2;
        X = reshape(pts(:,1), size(u));
        Y = reshape(pts(:,2), size(u));
        Z = reshape(pts(:,3), size(u));
        h = surf(X, Y, Z, 'FaceColor', color, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    end

localAxes = zeros(numLayers, 3);
localDeltas = zeros(numLayers, 1);

% 第一组：层1-3，相同光轴
% 拟合A' for layers 1-3 using P0, P1, P2
A_prime1 = fit_plane_normal(poincarePoints(1,2:4), poincarePoints(2,2:4), poincarePoints(3,2:4));
localAxes(1:3, :) = repmat(A_prime1, 3, 1);

% 计算delta for layer 1: using P0, P1, A1
localDeltas(1) = calculate_local_retardation(poincarePoints(1,2:4), poincarePoints(2,2:4), A_prime1);

% delta for layer 2: P1, P2
localDeltas(2) = calculate_local_retardation(poincarePoints(2,2:4), poincarePoints(3,2:4), A_prime1);

% delta for layer 3: P2, P3
localDeltas(3) = calculate_local_retardation(poincarePoints(3,2:4), poincarePoints(4,2:4), A_prime1);

% 第二组：层4-6
% 拟合A' for layers 4-6 using P3, P4, P5 (assuming P6 is P_{numLayers+1})
A_prime2_raw = fit_plane_normal(poincarePoints(4,2:4), poincarePoints(5,2:4), poincarePoints(6,2:4));

% 旋转校正：A2 = rotate A_prime2_raw by -sum(delta1:3) around A1
total_delta_upper = sum(localDeltas(1:3));
A2 = rotate_vector(A_prime2_raw, A_prime1, -total_delta_upper);
localAxes(4:6, :) = repmat(A2, 3, 1);

% 计算delta for layer 4: P3, P4, A2
localDeltas(4) = calculate_local_retardation(poincarePoints(4,2:4), poincarePoints(5,2:4), A2);

% delta for layer 5: P4, P5
localDeltas(5) = calculate_local_retardation(poincarePoints(5,2:4), poincarePoints(6,2:4), A2);

% delta for layer 6: P5, P6
localDeltas(6) = calculate_local_retardation(poincarePoints(6,2:4), poincarePoints(7,2:4), A2);

%% 可视化局部光轴
for i = 1:numLayers
    % 绘制局部光轴（箭头）
    idx = i + 1; % P_i at poincarePoints(idx, :)
    arrowLen = 0.25; % 统一箭头长度
    la = localAxes(i,:);
    la = la / (norm(la) + eps);
    ah = quiver3(poincarePoints(idx, 2), poincarePoints(idx, 3), poincarePoints(idx, 4), ...
        la(1)*arrowLen, la(2)*arrowLen, la(3)*arrowLen, 0, 'Color',[0 0.45 0.74], 'LineWidth', 1.6, 'MaxHeadSize', 1);
    % 用不同透明度和边框突出箭尾（点）
    scatter3(poincarePoints(idx,2), poincarePoints(idx,3), poincarePoints(idx,4), 72, 'MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','k');
    set(ah,'AutoScale','off');
end

% 绘制拟合平面 (A_prime1 & A2) 以半透明网格显示
center1 = mean(poincarePoints(1:3, 2:4), 1);
h1 = plot_plane(center1, A_prime1, [0.85, 0.33, 0.1]);
center2 = mean(poincarePoints(4:7, 2:4), 1);
h2 = plot_plane(center2, A2, [0.2, 0.6, 0.2]);

% 标注拟合法向量 A_prime1 和 A2（用长箭头显示）
quiver3(center1(1), center1(2), center1(3), A_prime1(1)*0.6, A_prime1(2)*0.6, A_prime1(3)*0.6, 0, 'Color',[0.8 0.3 0.1], 'LineWidth', 2, 'MaxHeadSize', 1.2);
text(center1(1)+0.03, center1(2)+0.03, center1(3)+0.03, 'A''_{group1}','Color',[0.8 0.3 0.1],'FontWeight','bold');
quiver3(center2(1), center2(2), center2(3), A2(1)*0.6, A2(2)*0.6, A2(3)*0.6, 0, 'Color',[0.2 0.6 0.2], 'LineWidth', 2, 'MaxHeadSize', 1.2);
text(center2(1)+0.03, center2(2)+0.03, center2(3)+0.03, 'A_{group2}','Color',[0.2 0.6 0.2],'FontWeight','bold');

% 图例与说明（把 legend 放外侧，避免遮挡）
legend([trajH], {'Poincaré trajectory'}, 'Location','bestoutside');
% 添加说明文字
annotation('textbox',[0.03 0.85 0.25 0.12],'String',sprintf('Local deltas (rad):\n%s',mat2str(localDeltas,3)),'FitBoxToText','on','BackgroundColor',[1 1 1 0.8]);

% 显示局部相位延迟
disp('局部相位延迟（单位：弧度）：');
disp(localDeltas);
% 如果出现 NaN（数值退化），进行提示并用 0 替代以便后续可视化
if any(isnan(localDeltas))
    warning('发现 NaN 的局部相位延迟——可能是退化点（P 与 A".prime 共线）。已将 NaN 替换为 0 以便查看结果。');
    localDeltas(isnan(localDeltas)) = 0;
    disp('已替换后的局部相位延迟（单位：弧度）：');
    disp(localDeltas);
end

%% 保存图像
outputPath = fullfile(pwd, 'output'); % 输出路径
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end
saveas(gcf, fullfile(outputPath, 'fig1_corrected.png')); % 保存为PNG格式
% 保存数据以便测试/验证
save(fullfile(outputPath,'fig1_results.mat'),'localDeltas','poincarePoints','localAxes');