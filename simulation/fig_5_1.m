%%% gen QUV from AM and phase with adding noise
%%% set the SNR to add the noise
%%% figure 5-1: result will change due to stochastic processes
%
% 说明：
% 本脚本在 Poincaré 球（单位 QUV 球）上生成一组带噪声的 Stokes 向量（QUV），
% 模拟随层累积的小相位旋转，并通过局部窗口拟合计算局部光轴（LA）和相位延迟（phase retardation），
% 最后将结果可视化并保存为图片和 FIG 文件。
%
% 主要步骤概览：
% 1) 初始化绘图（单位球）
% 2) 定义背景光轴（LA）
% 3) 生成以 p1 为中心的带高斯角噪点云（调用 generatePointsOnSphereSurface）
% 4) 对每个点施加累积相位旋转以模拟层间相位变化
% 5) 对数据做平滑处理并用滑动窗口调用 cumLAphr 计算局部 LA 与相位
% 6) 将结果绘制并保存到 output 文件夹

%% define the start point in PoinCare sphere
close all; clear; clc;
fig = figure; % 创建绘图窗口句柄，后面用于保存

% 绘制单位球作为 Poincaré 球的可视化参考
[X, Y, Z] = sphere(20);
SS = surf(X, Y, Z);
SS.FaceAlpha = 0.05; % 球体透明度
SS.EdgeAlpha = 0.05; % 网格线透明度
xlim([-1.1 1.1])
axis equal
xlabel('Q')
ylabel('U')
zlabel('V')
grid on; hold on

% LA BG: 背景线性偏振轴（作为旋转轴参考）
a1 = [1 0 0];
a1 = a1 / sqrt(sum(a1.^2));


% slab 1 stokes: 设置每层的小相位累积值（弧度）
phr1 = pi/40; % 每层的相位增量，用于模拟累积相位

% 目标方向 p1（作为点云中心）并归一化
p1 = [0 1 0];
p1 = p1 / sqrt(sum(p1.^2));
Np = 100; % 生成点数

% 生成以 p1 为中心、角扰动标准差为 0.2 的球面点云
p_slab1_awgn = generatePointsOnSphereSurface(Np, p1, 0.2);

% 构造一个以原始点为起点，并逐点按背景轴 a1 施加累积旋转的点集合
% 注意：这里先把第一个点复制为初始化，然后在循环中追加旋转后的点
p_slab1_awgn_rotBG(1, :) = p_slab1_awgn(1, :);
for i = 1:Np
    inpp = p_slab1_awgn(i, :);
    % 使用 rotationVectorToMatrix 对向量做旋转：旋转向量为 i*phr1*a1
    outp = rotationVectorToMatrix(i * phr1 * a1) * inpp';
    p_slab1_awgn_rotBG = [p_slab1_awgn_rotBG; outp'];
end

% 绘制原始（带噪）点与经累积相位旋转后的点，便于比较
plot3(p_slab1_awgn(:, 1), p_slab1_awgn(:, 2), p_slab1_awgn(:, 3), 'or', 'markersize', 3)
plot3(p_slab1_awgn_rotBG(:, 1), p_slab1_awgn_rotBG(:, 2), p_slab1_awgn_rotBG(:, 3), 'sb', 'markersize', 3)
%% fit LA (局部拟合光轴并计算相位延迟)
avgNum = 5; % 滑动窗口大小：用于局部拟合

% 对原始点和旋转后点分别做平滑处理以降低噪声（高斯滤波窗口长度不同）
p_slab1_awgn_smooth = smoothdata(p_slab1_awgn, 1, 'gaussian', 10);
p_slab1_awgn_smooth1 = smoothdata(p_slab1_awgn, 1, 'gaussian', 10);
% 归一化平滑后的向量，确保仍在单位球上
p_slab1_awgn_smooth1 = p_slab1_awgn_smooth1 ./ sqrt(sum(p_slab1_awgn_smooth1.^2, 2));
plot3(p_slab1_awgn_smooth1(:, 1), p_slab1_awgn_smooth1(:, 2), p_slab1_awgn_smooth1(:, 3), 'sm', 'markersize', 3)

% 对旋转后的点序列进行更强的平滑处理（用于对比）
p_slab1_awgn_rotBG_smooth = smoothdata(p_slab1_awgn_rotBG, 1, 'gaussian', 80);
p_slab1_awgn_rotBG_smooth1 = smoothdata(p_slab1_awgn_rotBG, 1, 'gaussian', 30);
p_slab1_awgn_rotBG_smooth1 = p_slab1_awgn_rotBG_smooth1 ./ sqrt(sum(p_slab1_awgn_rotBG_smooth1.^2, 2));
plot3(p_slab1_awgn_rotBG_smooth1(:, 1), p_slab1_awgn_rotBG_smooth1(:, 2), p_slab1_awgn_rotBG_smooth1(:, 3), 'sk-', 'markersize', 3)
% 归一化备用平滑结果
p_slab1_awgn_rotBG_smooth = p_slab1_awgn_rotBG_smooth ./ sqrt(sum(p_slab1_awgn_rotBG_smooth.^2, 2));

% 使用滑动窗口对邻近点做局部平面拟合和相位计算
for i = 1:Np - avgNum
    % 对四类数据：原始、旋转后、原始平滑、旋转后平滑，分别计算 LA 和相位
    [LAcum_raw, phr_raw] = cumLAphr(p_slab1_awgn(i:i + avgNum, :));
    [LAcum_rot, phr_rot] = cumLAphr(p_slab1_awgn_rotBG(i:i + avgNum, :));
    [LAcum_raw_sm, phr_raw_sm] = cumLAphr(p_slab1_awgn_smooth(i:i + avgNum, :));
    [LAcum_rot_sm, phr_rot_sm] = cumLAphr(p_slab1_awgn_rotBG_smooth(i:i + avgNum, :));

    % 绘制局部 LA 方向（作为短线段）以便可视化光轴方向
    if LAcum_raw_sm(2) > 0
        plot3([-LAcum_raw_sm(1), LAcum_raw_sm(1)], [-LAcum_raw_sm(2), LAcum_raw_sm(2)], [-LAcum_raw_sm(3), LAcum_raw_sm(3)], 'r-', 'linewidth', 0.1)
    else
        plot3([-LAcum_raw_sm(1), LAcum_raw_sm(1)], [-LAcum_raw_sm(2), LAcum_raw_sm(2)], [-LAcum_raw_sm(3), LAcum_raw_sm(3)], 'r-', 'linewidth', 0.1)
    end
    plot3([-LAcum_rot_sm(1), LAcum_rot_sm(1)], [-LAcum_rot_sm(2), LAcum_rot_sm(2)], [-LAcum_rot_sm(3), LAcum_rot_sm(3)], 'b-', 'linewidth', 0.1)

    % 将每个窗口结果保存到数组，便于后续统计或导出
    LAcums_raw(i, :) = LAcum_raw;
    LAcums_rot(i, :) = LAcum_rot;
    LAcums_raw_sm(i, :) = LAcum_raw_sm;
    LAcums_rot_sm(i, :) = LAcum_rot_sm;
    phrs_raw(i, :) = phr_raw;
    phrs_rot(i, :) = phr_rot;
    phrs_raw_sm(i, :) = phr_raw_sm;
    phrs_rot_sm(i, :) = phr_rot_sm;
end


%% 



function [LAcum, phR] = cumLAphr(planeData)
% cumLAphr - 在一个局部窗口的 QUV 点上拟合平面并计算局部光轴（LA）与相位延迟
%
% 输入:
%   planeData : M x 3 矩阵，M = Avnum+1，表示连续层的 QUV 向量（每行是 [Q U V]）
% 输出:
%   LAcum : 1x3 局部光轴（平面法向量，已归一化）
%   phR : 标量，相邻层之间的相位延迟估计（弧度）
%
% 方法概述：
% 1) 若数据包含非法值（如 Q 等于 0）则直接返回 0，以跳过该窗口
% 2) 对窗口中的点做去均值化后进行 SVD，第三列特征向量对应拟合平面的法向量
% 3) 对法向量的符号进行检查并修正（因为 SVD 得到的特征向量方向不唯一）
% 4) 计算每对相邻层的切向量（cross(Sin, LA)）并用它们的夹角的中位数来估计总体相位延迟

Avnum = size(planeData, 1) - 1; % 注意 planeData 行数 = Avnum+1

% 把 Q U V 分开便于处理
x1 = squeeze(planeData(:, 1));
y1 = squeeze(planeData(:, 2));
z1 = squeeze(planeData(:, 3));

% 检查是否存在 0 值（这里的 x1==0 意味着 Q 分量全部为 0 或存在 0）
nx = sum(sum(x1 == 0));
if x1 == 0 | nx > 0 % 如果 Q 分量有问题则跳过计算
    LAcum(1:3) = 0;
    phR = 0;
    return
end

% 去除窗口均值：把点云中心化以便进行 SVD 平面拟合
centeredPlane = planeData - mean(planeData);

% 用 SVD 对中心化点做拟合，第三列特征向量（V1(:,3)）为最小奇异值对应的方向，即平面法向量
[~, ~, V1] = svd(centeredPlane);
a1 = V1(1, 3);
b1 = V1(2, 3);
c1 = V1(3, 3);

% 归一化法向量分量，防止数值问题
normFactor = sqrt(a1.^2 + b1.^2 + c1.^2);
if normFactor ~= 0
    a1 = a1 ./ normFactor;
    b1 = b1 ./ normFactor;
    c1 = c1 ./ normFactor;
else
    a1 = 0; b1 = 0; c1 = 0;
end
% 防止 NaN（当 normFactor 为 0 或非常小时）
a1(isnan(a1)) = 0; b1(isnan(b1)) = 0; c1(isnan(c1)) = 0;

% check the axis orientation and reverse **十分重要
% SVD 得到的特征向量方向不确定（可正可负），需要用窗口首尾点判断方向是否反向
Sin1 = [x1(1), y1(1), z1(1)];
if x1(1) == 0
    % 如果第一个点的 Q 分量为 0，则使用第二个点作为起点参考
    Sin1 = [x1(2), y1(2), z1(2)];
end
Sout1 = [x1(Avnum), y1(Avnum), z1(Avnum)];
ax = [a1, b1, c1];

% 计算 Sin1 在 ax 上的投影并得到垂直分量，用于构造两个切向量
dd1 = dot(Sin1, ax);
SS1 = Sin1 - dd1 * ax; % Sin1 在平面上的分量

dd2 = dot(Sout1, ax);
SS2 = Sout1 - dd2 * ax; % Sout1 在平面上的分量

% 通过 cross(SS1, SS2) 与 ax 的点积判断方向
dd3 = cross(SS1, SS2);
dd4 = dot(dd3, ax);
if dd4 < 0
    % 如果方向不一致，需翻转法向量符号，使方向与点序列的“手性”一致
    a1 = -a1; b1 = -b1; c1 = -c1;
end
LAcum = [a1, b1, c1]; % 返回局部光轴（单位向量）

% cal phase retardation：对每对相邻层计算切向量之间的夹角余弦
val3 = zeros(1, Avnum); % 预分配以防未定义索引
for k = 1:Avnum
    Sin = [x1(k), y1(k), z1(k)];
    A = [a1, b1, c1];
    Sout = [x1(k + 1), y1(k + 1), z1(k + 1)];
    val1 = cross(Sin, A); % 切向量（Sin 在平面上的切向量）
    lval1 = sqrt(sum(val1 .^ 2));
    val2 = cross(Sout, A);
    lval2 = sqrt(sum(val2 .^ 2));
    if (lval1 == 0) || (lval2 == 0)
        val3(k) = 0; % 若任一切向量为 0（点与光轴共线），则无法计算夹角，置 0
    else
        % 两个切向量的点积除模长得到夹角余弦
        val3(k) = dot(val1, val2) / lval1 / lval2;
    end
    % 数值稳定：将余弦限制在 [-1, 1]
    val3(k) = max([val3(k), -1]);
    val3(k) = min([val3(k), 1]);
end

% 采用中位数以抑制离群值的影响，然后取 acos 得到相位延迟
phR = acos(median(val3));
end


function points_rot = generatePointsOnSphereSurface(numPoints, A, B)
% generatePointsOnSphereSurface - 在单位球面上围绕方向 A 生成带高斯角噪的点
%
% 输入：
%   numPoints - 生成点数（可选，默认 100）
%   A - 目标方向（三元向量），点云将以该方向为中心分布
%   B - 噪声标准差（以弧度为单位，用于生成小角度扰动）
% 输出：
%   points_rot - numPoints x 3 矩阵，返回在单位球面上的点坐标（每行是一个 [Q U V] 向量）
%
% 实现思路（两步变换）：
% 1) 先在北极方向 A0=[0 0 1] 周围生成一组小角度扰动点（扰动轴由随机 theta 决定）
% 2) 把这些扰动点整体绕某旋转轴 LA 旋转 Theta，使其基准方向从 A0 映射到目标方向 A

if nargin == 0
    numPoints = 100;
    A = [1 1 0];
    B = 0.2;
end

% 归一化目标方向，保证为单位向量
A = A / sqrt(sum(A .^ 2));

% 角偏移：每个点一个小的高斯随机角度（弧度）
disRad = B * randn(numPoints, 1); % normal distribution of small-angle perturbations
theta = 2 * pi * rand(numPoints, 1); % 环向角，用于生成不同方向的扰动轴

% 基准方向（北极），我们先在其附近生成扰动点再整体映射到 A
A0 = [0 0 1];
points = zeros(numPoints, 3); % 预分配，避免循环内动态扩展
for i = 1:numPoints
    % 构造一个与 A0 不同的随机参考向量 A1（保证不会与 A0 平行）
    A1 = [sin(theta(i)), cos(theta(i)), 1];
    orth = cross(A0, A1);
    % 若 orth 为 0 则跳过归一化（理论上 A1 不应与 A0 共线，但增加健壮性）
    if norm(orth) == 0
        orth = [1 0 0];
    else
        orth = orth / norm(orth);
    end
    % 使用 rotationVectorToMatrix 将 z 轴点小幅度旋转到扰动方向
    points(i, :) = (rotationVectorToMatrix(disRad(i) * orth) * [0 0 1]')';
end

% 计算把 A0 旋转到目标 A 所需的旋转角 Theta 与旋转轴 LA
% 注意：当 A 与 A0 平行时，LA 可能为 0 向量，需处理
Theta = acos(dot(A, A0));
LA = cross(A, A0);
normLA = sqrt(sum(LA .^ 2));
if normLA == 0
    % A 与 A0 平行或反向：旋转轴不可定义，直接把 points 作为最终结果（或对反向做 180° 旋转）
    points_rot = points;
    if dot(A, A0) < 0
        % A 与 A0 反向：需要绕任一垂直轴旋转 pi（180度）
        % 这里选择 x 轴作为旋转轴
        R = rotationVectorToMatrix(pi * [1 0 0]);
        points_rot = (R * points')';
    end
    return
end
LA = LA / normLA; % 归一化旋转轴

% 将北极扰动点整体旋转到以 A 为中心的分布
points_rot = zeros(numPoints, 3); % 预分配输出
for j = 1:numPoints
    inp = points(j, :);
    outp = rotationVectorToMatrix(Theta * LA) * inp';
    points_rot(j, :) = outp';
end
end

% 增加输出路径控制
outputPath = fullfile(pwd, 'output'); % 默认输出路径为当前目录下的 output 文件夹
if ~exist(outputPath, 'dir')
    mkdir(outputPath); % 如果路径不存在，则创建
end

% 保存图形
saveas(fig, fullfile(outputPath, 'figure5_1.png')); % 保存为 PNG 格式
savefig(fig, fullfile(outputPath, 'figure5_1.fig')); % 保存为 MATLAB FIG 格式