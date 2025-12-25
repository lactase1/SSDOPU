function [LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw] = FreeSpace_PSOCT_Optimized(Qm, Um, Vm, Stru, test_seg_top, h1, h2, Avnum, dopuMap, enableDopuPhaseSupp)
% FreeSpace_PSOCT_Optimized (修正版)
% 修复了遗漏的低 DOPU 区域预平滑逻辑

    %% 参数默认值与初始化
    if nargin < 10, enableDopuPhaseSupp = true; end
    if nargin < 9 || isempty(dopuMap), dopuMap = []; end
    
    [nZ, nX] = size(Qm);
    
    % 常量定义
    stLayer = 1;
    DTheta = 1.0675 * 2 * pi;
    dTheta = DTheta / 320;
    LABG = [1 0 0];
    
    % DOPU 阈值自动计算
    dopuLowThresh = 0.5;
    dopuMinScale = 0.1;
    if ~isempty(dopuMap)
        validDopu = dopuMap(dopuMap > 0 & ~isnan(dopuMap));
        if ~isempty(validDopu)
            dopuLowThresh = min(max(median(validDopu) * 1.2, 0.1), 0.8);
        end
    end

    %% Step 1: 数据平展
    StruF = zeros(size(Stru));
    QmF = zeros(size(Qm)); 
    UmF = zeros(size(Um)); 
    VmF = zeros(size(Vm));
    % 预分配 DOPU 缓存为完整帧大小，避免循环中动态扩展
    DopuF = zeros(nZ, nX);
    if ~isempty(dopuMap)
        % 仅拷贝传入的有效区域（保持与原逻辑一致）
        DopuF(1:size(dopuMap,1), 1:size(dopuMap,2)) = dopuMap;
    end

    for j = 1:nX
        top = test_seg_top(j);
        validLen = nZ - top + 1;
        if validLen > 0
            idxSrc = top:nZ;
            idxDst = 1:validLen;
            StruF(idxDst, j) = Stru(idxSrc, j);
            QmF(idxDst, j)   = Qm(idxSrc, j);
            UmF(idxDst, j)   = Um(idxSrc, j);
            VmF(idxDst, j)   = Vm(idxSrc, j);
            % 直接从已预分配 DopuF 读取，无需再修改尺寸
            % 如果没有输入 map，该处为 0，会在后续被忽略或触发回退计算
        end
    end

    %% Step 2: Stokes 矢量归一化与平滑
    normFactor = sqrt(QmF.^2 + UmF.^2 + VmF.^2);
    normFactor(normFactor == 0) = 1;
    
    SmapF = cat(3, QmF./normFactor, UmF./normFactor, VmF./normFactor);
    SmapF(isnan(SmapF)) = 0;
    SmapF1 = smoothdata(SmapF, 1, 'movmean', 16);

    %% Step 3: DDG 核心算法
    output_depth = nZ - 20; % 固定深度为总深度-20，避免尺寸不同
    cumLA_bg = zeros(output_depth, nX, 3);
    phR = zeros(output_depth, nX);
    phR_rmBG = zeros(output_depth, nX);
    phR_raw = zeros(output_depth, nX);
    cumLA_raw = zeros(output_depth, nX, 3);
    
    % 记录窗口DOPU用于后续平滑
    phd = zeros(output_depth, nX); 
    
    Stokes_DeDC = zeros(Avnum, 3);
    Stokes_DeMod = zeros(Avnum, 3);
    
    for j = 1:nX 
        surface_idx = test_seg_top(j);
        
        for i = 1:output_depth
            windowData = squeeze(SmapF1(i:i+Avnum, j, :)); 
            
            if any(windowData(:,1) == 0)
                continue;
            end
            
            % 优化 DOPU 使用逻辑：优先使用传入的 dopuMap（若存在且有效），只有在缺失或非法时才计算窗口 DOPU
            if ~isempty(dopuMap) && i <= size(dopuMap, 1) && j <= size(dopuMap, 2)
                winDOPU = dopuMap(i, j);
                % 若传入值为 NaN 或非法（<0）则回退到窗口估计
                if isnan(winDOPU) || winDOPU < 0
                    winDOPU = calcWindowDOPU(windowData);
                end
            else
                % 没有外部 DOPU map 时在本地计算
                winDOPU = calcWindowDOPU(windowData);
            end
            % 保证数值稳定且在 [0,1] 范围内
            if isnan(winDOPU), winDOPU = 0; end
            winDOPU = min(max(winDOPU, 0), 1);
            phd(i, j) = winDOPU; % 记录下来

            % 旋转补偿
            base_angle_idx = i + surface_idx + stLayer;
            for k = 1:Avnum
                vec_in = windowData(k, :);
                angle_dc = -1 * base_angle_idx * dTheta;
                Stokes_DeDC(k, :) = (rotationVectorToMatrix(angle_dc * LABG) * vec_in')';
                angle_mod = -1 * (base_angle_idx + k - 2) * dTheta;
                Stokes_DeMod(k, :) = (rotationVectorToMatrix(angle_mod * LABG) * vec_in')';
            end
            
            % --- 采用自适应多点拟合并使用切面圆约束（planeFitConstrained） ---
            % 动态确定拟合点数：DOPU 越高使用越少点，DOPU 越低使用越多点以增强稳定性
            opts.minFit = max(3, min(6, Avnum)); % 最小拟合点数（至少3）
            opts.maxFit = Avnum;                 % 最大拟合点数
            opts.r_min = 0.5;                    % 切面圆最小半径阈值（可调）

            currDOPU = min(max(winDOPU, 0), 1);
            nFit = round(opts.maxFit - (opts.maxFit - opts.minFit) * currDOPU);
            nFit = max(opts.minFit, min(opts.maxFit, nFit));

            % 取窗口头部的 nFit 个点用于拟合（保留深度方向信息）
            S_use_mod = Stokes_DeMod(1:nFit, :);
            S_use_dc  = Stokes_DeDC(1:nFit, :);

            % 对调制信号进行切面圆约束拟合
            [B_raw_vec, c_raw, r_raw] = planeFitConstrained(S_use_mod, opts.r_min);
            % 对去直流信号进行拟合（用于补偿/备选）
            [B_dc_vec, c_dc, r_dc] = planeFitConstrained(S_use_dc, opts.r_min);

            % [B_raw_vec, ~] = threePointFitting(S_use_mod);
            % [B_dc_vec, ~]  = threePointFitting(S_use_dc);

            % 方向符号校正：使用首中末点粗略叉乘确定旋转手性
            if size(S_use_mod,1) >= 3
                ref_dir = cross(S_use_mod(2,:) - S_use_mod(1,:), S_use_mod(end,:) - S_use_mod(round(end/2), :));
                if dot(B_raw_vec, ref_dir) < 0, B_raw_vec = -B_raw_vec; end
            end

            % 使用约束后的向量计算相位和原始累积 LA
            phR_raw(i, j) = calcPhase(S_use_mod, B_raw_vec);
            cumLA_raw(i, j, :) = -B_raw_vec;

            % 基本的窗口跨度度量（使用拟合的子集）
            P1 = S_use_mod(1, :);
            P3 = S_use_mod(end, :);
            T3_nor = norm(P3 - P1);

            % B_dedc 用于后续 DeltaB 逻辑，使用约束后的去直流结果
            B_dedc = B_dc_vec;
            
            if T3_nor < 0.065
                final_B = B_dedc;
                if final_B(1) > 0, final_B = -final_B; end 
                
                ph_temp = calcPhase(Stokes_DeDC, final_B);
                if i > 1 && phR(i-1, j) > 0
                    ph_temp = (ph_temp + phR(i-1, j)) / 2;
                end
                
                phR(i, j) = ph_temp;
                phR_rmBG(i, j) = min(max(ph_temp - dTheta, 0), 0.2);
            else
                % 复现 deltaB 逻辑
                mid = round((1+Avnum)/2);
                T1p = Stokes_DeDC(mid,:) - Stokes_DeDC(1,:);
                T2p = Stokes_DeDC(end,:) - Stokes_DeDC(mid,:);
                Tb1 = P1 - Stokes_DeDC(1,:);
                Tb2 = P3 - Stokes_DeDC(end,:);
                
                deltaB = cross(T1p, Tb2) + cross(T2p, Tb1) - cross(Tb1, Tb2);
                final_B = B_dedc + deltaB;
                final_B = final_B / norm(final_B);
                
                phR(i, j) = calcPhase(Stokes_DeMod, final_B);
                phR_rmBG(i, j) = phR(i, j);
            end
            
            % DOPU 抑制
            if enableDopuPhaseSupp && winDOPU < dopuLowThresh
                scale = max((winDOPU / dopuLowThresh)^2.0, dopuMinScale);
                phR_raw(i, j)  = phR_raw(i, j) * scale;
                phR(i, j)      = phR(i, j) * scale;
                phR_rmBG(i, j) = min(phR_rmBG(i, j), phR(i, j));
            end
            cumLA_bg(i, j, :) = -final_B;
        end
    end

    %% Step 3.5: [补回遗漏逻辑] 低 DOPU 区域预平滑
    % 这是之前遗漏的关键步骤，它使用小高斯核对低信噪比区域进行预处理
    if enableDopuPhaseSupp
        % 定义小高斯核 (与原版一致)
        smallGauss = fspecial('gaussian', [5 5], 1.2);
        
        % 对三种相位图进行滤波
        phR_s      = imfilter(phR, smallGauss, 'replicate');
        phR_raw_s  = imfilter(phR_raw, smallGauss, 'replicate');
        phR_rmBG_s = imfilter(phR_rmBG, smallGauss, 'replicate');
        
        % 构建掩码：DOPU 有效且小于阈值
        maskLow = (phd > 0) & (phd < dopuLowThresh);
        
        % 仅替换低 DOPU 区域
        if any(maskLow(:))
            phR(maskLow)      = phR_s(maskLow);
            phR_raw(maskLow)  = phR_raw_s(maskLow);
            phR_rmBG(maskLow) = phR_rmBG_s(maskLow);
        end
    end

    %% Step 4: 空间滤波
    [cumLA_bg_gF] = filterVectorField(cumLA_bg, h2);
    [cumLA_raw_gF] = filterVectorField(cumLA_raw, h2);
    
    phR_gF      = imfilter(phR, h2, 'replicate');
    phR_rmBG_gF = imfilter(phR_rmBG, h2, 'replicate');
    phR_raw_gF  = imfilter(phR_raw, h2, 'replicate');

    %% Step 5: 三维光轴递推
    bfloaxis3D     = drLA_Optimized(cumLA_bg_gF, -phR_gF, h2);
    bfloaxis3D_raw = drLA_Optimized(cumLA_raw_gF, -phR_raw_gF, h2);

    %% Step 6: 结果还原
    [LA, PhR, cumLA] = unflatten(bfloaxis3D, phR_rmBG_gF, cumLA_bg_gF, test_seg_top, nZ);
    [LA_raw, PhR_raw, cumLA_raw] = unflatten(bfloaxis3D_raw, phR_raw_gF, cumLA_raw_gF, test_seg_top, nZ);
end

%% ================= 子函数区域 =================
% (保持不变，与上一版一致)
function dopu = calcWindowDOPU(data)
    avgS = mean(data, 1);
    dopu = sqrt(sum(avgS.^2));
end

%% 多点 SVD 平面拟合（替代三点法，增加鲁棒性）
function [B, phR_median] = multiPointFitting(S)
    % 如果点数过少则退回三点法
    if size(S,1) < 3
        [B, phR_median] = threePointFitting(S);
        return;
    end

    % 去中心化
    center = mean(S, 1);
    Sc = S - center;

    % SVD 分解，最小奇异值对应的右奇异向量为平面法向量
    [~, ~, V] = svd(Sc, 0);
    B = V(:, end)';

    % 归一化与容错
    bn = norm(B);
    if bn < 1e-9
        B = [1 0 0];
    else
        B = B / bn;
    end

    % 方向符号校正：利用首中末三点的叉积作为参考方向
    P_start = S(1, :);
    P_mid = S(round(end/2), :);
    P_end = S(end, :);
    B_ref = cross(P_mid - P_start, P_end - P_mid);
    if dot(B, B_ref) < 0
        B = -B;
    end

    % 计算相位延迟（复用现有 calcPhase）
    phR_median = calcPhase(S, B);
end

function [B, phR_median] = threePointFitting(S)
    Avnum = size(S, 1);
    mid = round((1+Avnum)/2);
    P1 = S(1, :);
    P2 = S(mid, :);
    P3 = S(end, :);
    T1 = P2 - P1;
    T2 = P3 - P2;
    B = cross(T1, T2);
    bnorm = norm(B);
    if bnorm > 1e-6, B = B / bnorm; else, B = [1 0 0]; end
    phR_median = calcPhase(S, B);
end

function phR_val = calcPhase(S, B)
    Avnum = size(S, 1);
    vals = zeros(1, max(0, Avnum-2));
    for k = 1:max(0, Avnum-2)
        T1 = S(k+1,:) - S(k,:);
        T2 = S(k+2,:) - S(k+1,:);
        N1 = cross(B, T1);
        N2 = cross(B, T2);
        denom = norm(N1) * norm(N2);
        if denom > 1e-9
            val = dot(N1, N2) / denom;
            vals(k) = max(min(val, 1), -1);
        else
            vals(k) = 1;
        end
    end
    if isempty(vals)
        phR_val = 0;
    else
        phR_val = acos(median(vals));
    end
end

%% 切面圆约束实现（planeFitConstrained）
function [n, c, r] = planeFitConstrained(P, r_min)
    % P: [N x 3] 点云
    N = size(P,1);
    if N < 3
        % 退回到三点法（或简单 SVD）以保证稳定
        if N == 2
            v = cross(P(2,:)-P(1,:), [1 0 0]); vnorm = norm(v);
            if vnorm < 1e-6, n = [1 0 0]; else n = (v / vnorm); end
            c = dot(n, mean(P,1)');
            r = sqrt(max(0, 1 - c^2));
            return;
        else
            n = [1 0 0]; c = 0; r = 1; return;
        end
    end

    % 1) SVD 无约束拟合
    Pm = mean(P, 1);
    P_centered = P - Pm;
    [~, ~, V] = svd(P_centered, 0);
    n = V(:, 3)';
    if dot(n, Pm) < 0, n = -n; end
    % 2) 初始距离与半径
    c0 = dot(n, Pm);
    r0 = sqrt(max(0, 1 - c0^2));

    % 3) 约束检查与软调整
    if r0 >= r_min
        c = c0; r = r0; return;
    end

    % 软约束参数
    fit_error = sum((P * n' - c0).^2) / N;
    noise_level = fit_error / (1 + eps);
    soft_factor = 100;
    alpha = 1.0 / (1.0 + noise_level * soft_factor);
    alpha = max(0.3, min(0.9, alpha));

    c_max = sqrt(max(0, 1 - r_min^2));
    c = sign(c0) * (alpha * abs(c0) + (1 - alpha) * c_max);
    r = sqrt(max(0, 1 - c^2));
end

function [V_out] = filterVectorField(V_in, h)
    V1 = imfilter(V_in(:,:,1), h, 'replicate');
    V2 = imfilter(V_in(:,:,2), h, 'replicate');
    V3 = imfilter(V_in(:,:,3), h, 'replicate');
    normV = sqrt(V1.^2 + V2.^2 + V3.^2);
    normV(normV==0) = 1;
    V_out = cat(3, V1./normV, V2./normV, V3./normV);
    V_out(isnan(V_out)) = 0;
end

function [axiss2] = drLA_Optimized(loaxis22, phR_gF, ~)
    [nZ, nX, ~] = size(loaxis22);
    axiss2 = zeros(size(loaxis22));
    for j = 1:nX
        current_axis = squeeze(loaxis22(1, j, :))';
        current_rotM = rotationVectorToMatrix(-phR_gF(1,j)/2 * current_axis);
        axiss2(1, j, :) = current_axis;
        for i = 2:nZ
            a = squeeze(loaxis22(i, j, :));
            new_axis = (current_rotM * a)';
            new_axis = new_axis / (norm(new_axis) + 1e-9);
            axiss2(i, j, :) = new_axis;
            d1 = -phR_gF(i, j) / 2.0 * new_axis;
            current_rotM = rotationVectorToMatrix(d1) * current_rotM;
        end
    end
end

function [LA, PhR, cumLA] = unflatten(srcLA, srcPhR, srcCumLA, tops, nZ)
    [outZ, nX, ~] = size(srcCumLA);
    LA = zeros(outZ, nX, 3);
    PhR = zeros(outZ, nX);
    cumLA = zeros(outZ, nX, 3);
    for j = 1:nX
        t = tops(j);
        len = min(outZ, nZ - t + 1);
        if len > 0
            idx_src = 1:len;
            idx_dst = t:(t+len-1);
            if idx_dst(end) > size(LA, 1), idx_dst = idx_dst(idx_dst <= size(LA,1)); end
            len_act = length(idx_dst);
            LA(idx_dst, j, :) = srcLA(1:len_act, j, :);
            PhR(idx_dst, j) = srcPhR(1:len_act, j);
            cumLA(idx_dst, j, :) = srcCumLA(1:len_act, j, :);
        end
    end
end