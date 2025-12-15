% 本脚本用于PS-OCT后巩膜检测的核心算法实现，参考《Robust Detection of Posterior Scleral using PS-OCT》理论。
% 主要流程包括：QUV归一化、DDG法光轴与相位拟合、光轴与相位滤波、结果还原。
% 关键概念：
%   QUV：Stokes矢量的三个分量，反映偏振态信息。
%   DDG（Depth-Dependent Gradient）：通过拟合局部Stokes矢量变化，估算光轴方向和相位延迟。
%   光轴（LA）：局部双折射主轴方向，反映组织微结构。
%   相位延迟（PhR）：光在双折射介质中传播的相位差，反映组织厚度和性质。
%   平展处理：将A-scan数据按表面索引对齐，便于后续窗口处理。
%   背景去除：通过旋转矩阵消除系统性偏振背景。


function [LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw, LA_forced, PhR_forced] = FreeSpace_PSOCT_3_DDG_rmBG_7(Qm, Um, Vm, Stru, test_seg_top, h1, h2, Avnum, dopuMap, enableDopuPhaseSupp)
    if nargin < 10
        enableDopuPhaseSupp = true;
    end
    if nargin < 9 || isempty(dopuMap)
        dopuMap = [];
    end
    %% Step 1. 数据平展（A-scan表面对齐）
    % 目的：将每条A-scan的表面对齐，消除曲面影响，便于窗口化处理。
    % 参考论文中“曲面展开”思想。

    StruF = Stru * 0;
    QmF = Qm * 0;
    UmF = Um * 0;
    VmF = Vm * 0;
    if ~isempty(dopuMap)
        DopuF = Stru * 0;
    else
        DopuF = [];
    end

    for j = 1:size(Stru, 2)
        StruF(1:size(Stru, 1) - test_seg_top(j) + 1, j) = Stru(test_seg_top(j):size(Stru, 1), j);
        QmF(1:size(Stru, 1) - test_seg_top(j) + 1, j) = Qm(test_seg_top(j):size(Stru, 1), j);
        UmF(1:size(Stru, 1) - test_seg_top(j) + 1, j) = Um(test_seg_top(j):size(Stru, 1), j);
        VmF(1:size(Stru, 1) - test_seg_top(j) + 1, j) = Vm(test_seg_top(j):size(Stru, 1), j);
        if ~isempty(DopuF)
            maxRow = min(size(dopuMap, 1), size(Stru, 1));
            if test_seg_top(j) <= maxRow
                srcInd = test_seg_top(j):maxRow;
                destLen = numel(srcInd);
                DopuF(1:destLen, j) = dopuMap(srcInd, j);
            end
        end
    end

    %% Step 2. QUV归一化（Stokes矢量单位化）
    % 目的：将每个体素的Stokes矢量归一化到单位球面，便于后续DDG拟合。
    % 参考论文“Stokes矢量归一化”部分。

    QF = QmF;
    UF = UmF;
    VF = VmF;

    QF1 = QF ./ sqrt(QF.^2 + UF.^2 + VF.^2);
    UF1 = UF ./ sqrt(QF.^2 + UF.^2 + VF.^2);
    VF1 = VF ./ sqrt(QF.^2 + UF.^2 + VF.^2);
    
    %把Nan赋值成0
    QF1(isnan(QF1)) = 0;
    UF1(isnan(UF1)) = 0;
    VF1(isnan(VF1)) = 0;

    SmapF = cat(3, QF1, UF1, VF1); % 三维矩阵，每个体素为单位Stokes矢量
    SmapF1 = smoothdata(SmapF, 1, 'movmean', 16); % 深度方向平滑，抑制噪声

    % 强制半径法光轴与相位（DDG平行分支，统一返回供对比）
    forcedParam.r_min = 0.9; % 最小切圆半径约束（靠近1约束更强）
    forcedParam.output_depth_forced = max(size(SmapF1, 1) - Avnum, 1); % 与DDG输出深度对齐
    forcedParam.window_forced_max = Avnum + 1; % 复用当前窗口长度
    [LA_forced_flat, PhR_forced_flat] = computeFromQUV_Forced(SmapF1, forcedParam);

    %% Step 3. DDG法光轴与相位拟合
    % 目的：在Avnum窗口内，利用Stokes矢量变化，估算局部双折射主轴（光轴）和相位延迟。
    % 详细流程：
    % 1. 遍历每个B-scan的每个A-scan（横向j，纵向i），对每个深度窗口进行DDG拟合。
    % 2. 提取窗口内的单位Stokes矢量（planeData），判断数据有效性（是否有全零或无效点）。
    % 3. 对有效数据，分别进行背景去除（ps_bg_rdu）和深度补偿（ps_bg_rm）旋转，得到两个Stokes矢量序列。
    % 4. 用三点法拟合主轴：取窗口顶部、中部、底部三个点，计算两段向量T1、T2，叉积得主轴方向B_bg_rm。
    % 5. 用DDG切向法（ddgTsPhr）计算局部相位延迟（phR_raw），主轴方向（cumLA_raw）。
    % 6. 再用未去背景的三点法和补偿，得到Bp和deltaB，合成最终主轴B。
    % 7. 根据窗口跨度T3_nor，判断是否用去背景结果还是补偿结果计算相位。
    % 8. 最终将主轴和相位结果存入cumLA_bg和phR_rmBG。

    stLayer = 1;                % DDG拟合时的深度起始层
    DTheta = 1.0675 * 2 * pi;   % 总相位延迟范围（与仪器参数相关）
    dTheta = DTheta / 320;      % 每层的相位步进（与采样点数相关）
    LABG = [1 0 0];             % 参考主轴方向（背景主轴）
    phd = zeros(size(SmapF1, 1), size(SmapF1, 2)); % 用于记录窗口内的DOPU估计
    dopuLowThresh = 0.5; % DOPU阈值，低于该值触发相位下压（初始值，后续可调整）
    dopuMinScale = 0.1;   % 避免极端情况下完全归零，保留最小比例

    % 自动选择 DOPU 阈值（如果传入了 DOPU map）
    if ~isempty(dopuMap)
        dopuVals = dopuMap(:);
        dopuVals = dopuVals(~isnan(dopuVals));
        dopuVals = dopuVals(dopuVals > 0);
        if ~isempty(dopuVals)
            % 使用中位数的比例作为阈值，避免被极端最小值干扰
            medDopu = median(dopuVals);
            autoRatio = 1.2; % 可调整的比例系数
            autoThresh = medDopu * autoRatio;
            % 限定阈值到一个合理区间，防止过小或过大
            autoThresh = min(max(autoThresh, 0.1), 0.8);
            dopuLowThresh = autoThresh;
            fprintf('[DOPU自动阈值] median=%.3f -> dopuLowThresh=%.3f (ratio=%.2f)\n', medDopu, dopuLowThresh, autoRatio);
        end
    end

    for j = 1:size(SmapF1, 2)   % 遍历每个A-scan（横向）
        dirCheck = 1;           % 主轴方向修正标志
        for i = 1:size(SmapF1, 1) - Avnum % 遍历每个深度窗口
            planeData = squeeze(SmapF1(i:i+Avnum, j, :)); % 当前窗口内的单位Stokes矢量，大小[Avnum+1,3]
            x1 = squeeze(planeData(:, 1)); % 当前窗口所有点的Q分量
            nx = sum(x1 == 0);             % 统计Q分量为0的点数

            if any(x1 == 0) || nx > 0           % 判断窗口内是否有无效数据
                cumLA_bg(i, j, 1:3) = 0;       
                phR(i, j) = 0; 
                phR_rmBG(i, j) = 0; 
                phR_raw(i, j) = 0; 
                phd(i, j) = 0; 
                rotationVector(i, j, 1:3) = 0; 
            else
                if ~isempty(DopuF)
                    windowDOPU = DopuF(i, j);
                    if windowDOPU == 0
                        avgStokes = mean(planeData, 1); % 回退到窗口平均估计
                        windowDOPU = sqrt(sum(avgStokes .^ 2));
                    end
                else
                    avgStokes = mean(planeData, 1); % 用窗口平均Stokes估算DOPU
                    windowDOPU = sqrt(sum(avgStokes .^ 2));
                end
                windowDOPU = min(max(windowDOPU, 0), 1); % 限制到[0,1]
                phd(i, j) = windowDOPU; % 记录DOPU分布，便于后续分析
                ps_bg_rdu = []; % 存储去背景后的Stokes矢量序列
                ps_bg_rm = []; % 存储深度补偿后的Stokes矢量序列
                for jL = 1:Avnum
                    inpp = planeData(jL, :); % 当前窗口第jL层的Stokes矢量
                    outp = rotationVectorToMatrix(-1 * (i + test_seg_top(j) + stLayer - 0) * dTheta * LABG) * inpp'; % 去除背景后的旋转结果
                    outp2 = rotationVectorToMatrix(-1 * (i + test_seg_top(j) + stLayer + jL - 2) * dTheta * LABG) * inpp'; % 去除背景+深度补偿后的旋转结果
                    ps_bg_rdu = [ps_bg_rdu; outp']; % 累加去背景结果
                    ps_bg_rm = [ps_bg_rm; outp2']; % 累加深度补偿结果
                end

                P1 = ps_bg_rm(1, :);                      % 窗口顶部点的Stokes矢量
                P2 = ps_bg_rm(round((1 + Avnum) / 2), :); % 窗口中部点的Stokes矢量
                P3 = ps_bg_rm(Avnum, :);                  % 窗口底部点的Stokes矢量
                % ---【窗口三点法】---
                % 目的：用窗口顶部、中部、底部三个点，拟合当前局部组织的双折射主轴方向
                T1 = P2 - P1;                               % 顶部到中部的向量，反映Stokes矢量变化趋势
                T2 = P3 - P2;                               % 中部到底部的向量，反映Stokes矢量变化趋势
                B_bg_rm = cross(T1, T2);                    % 用两段变化向量的叉积，得到切平面的法向量，即主轴方向
                B_bg_rm = B_bg_rm / sqrt(sum(B_bg_rm .^ 2)); % 单位化主轴方向，保证为单位向量
                [Subphr_bg_rm_Ts] = ddgTsPhr(ps_bg_rm, B_bg_rm); % 用DDG切向法计算窗口内的相位变化，反映局部双折射特性
                phR_raw(i, j) = acos(median(Subphr_bg_rm_Ts)); % 取中位数作为窗口的相位延迟，抑制异常值影响
                cumLA_raw(i, j, 1:3) = -B_bg_rm; % 局部窗口的主轴方向，存入结果

                % ---【去背景三点法】---
                % 目的：用去除背景后的Stokes矢量，拟合主轴方向和补偿量，为后续主轴修正做准备
                P1p = ps_bg_rdu(1, :); % 去偏振调制后顶部点
                % P2p = ps_bg_rdu(2, :);
                % P3p = ps_bg_rdu(3, :); 
                P2p = ps_bg_rdu(round((1 + Avnum) / 2), :); 
                P3p = ps_bg_rdu(Avnum, :); 

                T1p = P2p - P1p; % 去偏振调制切向量1
                T2p = P3p - P2p; % 去偏振调制切向量2
                Tb1 = P1 - P1p; % 顶部点的补偿向量，反映背景影响
                Tb2 = P3 - P3p; % 底部点的补偿向量
                Bp = cross(T1p, T2p); % 未去背景主轴方向，用于后续主轴合成
                T3_nor = norm(P3 - P1); % 顶部到底部的距离（窗口跨度），用于判断是否需要补偿

                if T3_nor < 0.065 % 判断窗口跨度是否较小
                    deltaB = [0 0 0]; % 小跨度不做调制去除
                    if dirCheck && (Bp(1) > 0)
                        Bp = -Bp; % 修正主轴方向
                    end
                else
                    dirCheck = 0;
                    deltaB = cross(T1p, Tb2) + cross(T2p, Tb1) - cross(Tb1, Tb2); % 大跨度时补偿主轴
                end

                B = Bp + deltaB; % 合成最终主轴方向
                B = B / sqrt(sum(B .^ 2)); % 单位化主轴

                if T3_nor < 0.065 % 小跨度用去背景结果计算相位
                    phR(i, j) = acos(median(ddgTsPhr(ps_bg_rdu, B))); % 用去背景主轴计算相位延迟
                    if (i > 1) && (phR(i - 1, j) > 0)
                        phR(i, j) = (phR(i, j) + phR(i - 1, j)) / 2; % 与上一深度做一阶低通滤波，压制小窗口噪声
                    end
                    phR_rmBG(i, j) = max([phR(i, j) - 1 * dTheta, 0]); % 去除背景后的相位延迟
                    phR_rmBG(i, j) = min([phR(i, j), 0.2]); % 限制最大相位延迟
                else
                    phR(i, j) = acos(median(ddgTsPhr(ps_bg_rm, B))); % 大跨度用补偿主轴计算相位延迟
                    phR_rmBG(i, j) = phR(i, j);
                end

                if enableDopuPhaseSupp && (windowDOPU < dopuLowThresh) % DOPU较低表示偏振统一性差，抑制相位
                    % 强化衰减：使用幂次映射将低DOPU区域的缩放因子进一步减小
                    % 例如 gamma=2 会比线性缩放带来更强的抑制
                    gamma = 2.0;
                    scale = max((windowDOPU / dopuLowThresh) ^ gamma, dopuMinScale); % 计算衰减系数
                    phR_raw(i, j) = phR_raw(i, j) * scale; % 同步压制原始相位
                    phR(i, j) = phR(i, j) * scale; % 压制当前相位
                    phR_rmBG(i, j) = min(phR_rmBG(i, j), phR(i, j)); % 确保去背景相位不高于压制后数值
                end
                cumLA_bg(i, j, 1:3) = -B; % 存储最终主轴方向
            end
        end
    end

    %% 在应用大尺度结构先验滤波前，针对 DOPU 低的区域做局部高斯平滑（可减少伪相位）
    % 逻辑：先用 smallGaussian 对整个 phR/phR_raw/phR_rmBG 做一次平滑，
    % 然后只把 DOPU 低的像素替换为平滑后的值，保留高 DOPU 区域的原始精细结构。
    if exist('phd', 'var') && enableDopuPhaseSupp && exist('phR', 'var')
        % 小尺度高斯核（可调）
        smallGauss = fspecial('gaussian', [5 5], 1.2);
        phR_s = imfilter(phR, smallGauss, 'replicate');
        phR_raw_s = imfilter(phR_raw, smallGauss, 'replicate');
        phR_rmBG_s = imfilter(phR_rmBG, smallGauss, 'replicate');

        % 确保 phd 与 phR 大小一致（phd 原来为 size(SmapF1,1)，phR 为 size(SmapF1,1)-Avnum）
        rPhR = size(phR, 1);
        cPhR = size(phR, 2);
        phd_cropped = phd(1:rPhR, 1:cPhR);

        % 构建低 DOPU 掩码（phd 中 0 表示无效，忽略）
        maskLow = (phd_cropped > 0) & (phd_cropped < dopuLowThresh);
        nLow = nnz(maskLow);
        if nLow > 0
            fprintf('[DOPU平滑] 对 %d 个像素应用小尺度高斯平滑 (阈值=%.3f)\n', nLow, dopuLowThresh);
            phR(maskLow) = phR_s(maskLow);
            phR_raw(maskLow) = phR_raw_s(maskLow);
            phR_rmBG(maskLow) = phR_rmBG_s(maskLow);
        end
    end

    %% Step 4. 光轴与相位滤波（结构先验）
    % 目的：利用大尺度滤波（h2），结合组织结构先验，进一步抑制噪声，提升主轴与相位的空间一致性。

    [cumLA_bg_gF] = LAgFfilt(cumLA_bg, h2);
    [cumLA_raw_gF] = LAgFfilt(cumLA_raw, h2);
    phR_gF = imfilter(phR, h2, 'replicate');
    phR_rmBG_gF = imfilter(phR_rmBG, h2, 'replicate');
    phR_raw_gF = imfilter(phR_raw, h2, 'replicate');

    %% Step 5. 旋转矩阵与主轴递推
    % 目的：根据滤波后的主轴和相位，递推计算各深度的主轴方向，实现三维主轴场恢复。

    for j = 1:size(SmapF1, 2)
        ax_rot(1, :) = cumLA_bg_gF(1, j, :);
        rotationVector(1, j, :) = -phR_gF(1, j) / 2 * ax_rot;
        ax_rot(1, :) = cumLA_raw_gF(1, j, :);
        rotationVector_raw(1, j, :) = -phR_raw_gF(1, j) / 2 * ax_rot;
    end

    bfloaxis3D = drLA(cumLA_bg_gF, -phR_gF, rotationVector);
    bfloaxis3D_raw = drLA(cumLA_raw_gF, -phR_raw_gF, rotationVector_raw);

    %% Step 6. 曲面还原（结果映射回原始坐标）
    % 目的：将平展处理后的主轴和相位结果映射回原始曲面坐标，保证与原始数据一致。

    bfphrr = phR_rmBG_gF;
    bfloaxis3D(isnan(bfloaxis3D)) = 0;
    bfphrr(isnan(bfphrr)) = 0;
    bfphrr_raw = phR_raw_gF;
    bfloaxis3D_raw(isnan(bfloaxis3D_raw)) = 0;
    bfphrr_raw(isnan(bfphrr_raw)) = 0;

    LA = bfloaxis3D * 0;
    PhR = bfphrr * 0;
    cumLA = bfloaxis3D * 0;
    LA_raw = bfloaxis3D * 0;
    PhR_raw = bfphrr * 0;
    cumLA_raw = bfloaxis3D * 0;
    LA_forced = bfloaxis3D * 0;
    PhR_forced = bfphrr * 0;

    for j = 1:size(Stru, 2)
        toInd = test_seg_top(j):size(LA, 1);
        fromInd = 1:size(LA, 1) - test_seg_top(j) + 1;
        LA(toInd, j, :) = bfloaxis3D(fromInd, j, :);
        PhR(toInd, j) = bfphrr(fromInd, j);
        cumLA(toInd, j, :) = cumLA_bg_gF(fromInd, j, :);
        LA_raw(toInd, j, :) = bfloaxis3D_raw(fromInd, j, :);
        PhR_raw(toInd, j) = bfphrr_raw(fromInd, j);
        cumLA_raw(toInd, j, :) = cumLA_raw_gF(fromInd, j, :);
        % 强制半径法结果映射回原始曲面坐标
        if size(LA_forced_flat, 1) >= fromInd(end)
            LA_forced(toInd, j, :) = LA_forced_flat(fromInd, j, :);
        else
            LA_forced(toInd, j, :) = LA_forced_flat(end, j, :);
        end
        if size(PhR_forced_flat, 1) >= fromInd(end)
            PhR_forced(toInd, j) = PhR_forced_flat(fromInd, j);
        else
            PhR_forced(toInd, j) = PhR_forced_flat(end, j);
        end
    end
end

%% DDG切向法相位计算（参考论文公式）
function [Subphr_bg_rm_Ts] = ddgTsPhr(ps_bg_rm,B_bg_rm)
    Avnum = size(ps_bg_rm,1);
    for k = 1:Avnum-2
        T1_bg_rm = ps_bg_rm(k+1,:)-ps_bg_rm(k,:);
        T2_bg_rm = ps_bg_rm(k+2,:)-ps_bg_rm(k+1,:);
        N1_bg_rm = cross(B_bg_rm,T1_bg_rm); 
        N2_bg_rm = cross(B_bg_rm,T2_bg_rm);
        Subphr_bg_rm_Ts(k) = dot(N1_bg_rm,N2_bg_rm)/(sqrt(sum(N1_bg_rm.^2)*sum(N2_bg_rm.^2)));
    end
    Subphr_bg_rm_Ts(Subphr_bg_rm_Ts > 1) = 1; 
    Subphr_bg_rm_Ts(Subphr_bg_rm_Ts < -1) = -1;
end

%% DDG点向法相位计算（补充，部分场景使用）
function [Subphr_bg_rm_Ps] = ddgPsPhr(ps_bg_rm,B_bg_rm) %#ok<DEFNU>
    Avnum = size(ps_bg_rm,1);
    for k = 1:Avnum-1
        P1t_bg_rm = cross(ps_bg_rm(k,:),B_bg_rm);
        P2t_bg_rm = cross(ps_bg_rm(k+1,:),B_bg_rm);
        Subphr_bg_rm_Ps(k) = dot(P1t_bg_rm,P2t_bg_rm)/(sqrt(sum(P1t_bg_rm.^2)*sum(P2t_bg_rm.^2)));
    end
    Subphr_bg_rm_Ps(Subphr_bg_rm_Ps > 1) = 1; Subphr_bg_rm_Ps(Subphr_bg_rm_Ps < -1) = -1;
end

%% 光轴滤波（结构先验，空间一致性）
function [loaxis3] = LAgFfilt(loaxis2,h2)
    Temp_ax1(:,:) = loaxis2(:,:,1); 
    Temp_ax2(:,:) = loaxis2(:,:,2); 
    Temp_ax3(:,:) = loaxis2(:,:,3);
    Temp_ax1 = imfilter(Temp_ax1,h2,'replicate'); 
    Temp_ax2 = imfilter(Temp_ax2,h2,'replicate'); 
    Temp_ax3 = imfilter(Temp_ax3,h2,'replicate');
    Temp_ax1 = Temp_ax1./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
    Temp_ax2 = Temp_ax2./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
    Temp_ax3 = Temp_ax3./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
    Temp_ax1(isnan(Temp_ax1)) = 0; 
    Temp_ax2(isnan(Temp_ax2)) = 0; 
    Temp_ax3(isnan(Temp_ax3)) = 0;
    loaxis3 = cat(3,Temp_ax1,Temp_ax2,Temp_ax3);
end

%% 光轴递推（主轴场恢复，参考论文递推公式）
function [axiss2] = drLA(loaxis22,phR_gF,rotationVector2)
    for j = 1:size(loaxis22,2)
        for i = 1:size(loaxis22,1)
            vector1 = squeeze(rotationVector2(1,j,:));
            a = [loaxis22(i,j,1),loaxis22(i,j,2),loaxis22(i,j,3)];
            if i==1
                rotationMatrix1 = rotationVectorToMatrix(vector1);
                axis2(i,j,1:3) = a;
            else
                axis2(i,j,1:3) = rotationMatrix1*a';
                Qa = axis2(i,j,1); 
                Ua = axis2(i,j,2); 
                Va = axis2(i,j,3);
                d1 = -phR_gF(i,j)/2.0*[Qa,Ua,Va];
                rotationMatrix1 = rotationVectorToMatrix(d1)*rotationMatrix1;
                Q2 = Qa/sqrt(Qa^2+Ua^2+Va^2); 
                U2 = Ua/sqrt(Qa^2+Ua^2+Va^2); 
                V2 = Va/sqrt(Qa^2+Ua^2+Va^2);
                V2(isnan(V2)) = 0;
                axiss2(i,j,:) = cat(3,Q2,U2,V2);
            end
        end
    end
end
