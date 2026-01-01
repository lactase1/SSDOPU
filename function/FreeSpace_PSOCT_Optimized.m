function [LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw] = FreeSpace_PSOCT_Optimized(Qm, Um, Vm, Stru, test_seg_top, h1, h2, Avnum, dopuMap, enableDopuPhaseSupp, dopuThreshold, enableOutputAdaptive, kRL_output, kRU_output, outputDopuThreshold, enableBottomLayerPhaseReduction, bottomLayerDepth, adaptiveFilterBottomDepth, bottomDopuThreshold)
% FreeSpace_PSOCT_Optimized (双层融合版)
% 修复: 移除了导致垂直条纹的轴向惯性(Axial Inertia)逻辑，改用双层融合策略
% 主要改进:
%   1. 移除不稳定的 median(validDopu)*1.2 阈值自动计算
%   2. 直接使用预处理 DOPU Map，移除循环内冗余计算
%   3. 基于 DOPU 的自适应拟合策略（高 DOPU 小窗口，低 DOPU 大窗口+约束）
%   4. 支持输出端自适应滤波：根据DOPU值动态调整光轴/延迟的滤波强度
%   5. 支持底层相位延迟减小：在深层区域对低DOPU像素降低相位延迟（基于绝对深度与配置的阈值）
%   6. 支持分区域滤波：DOPU自适应滤波仅对底部指定层数生效 (基于绝对深度)
%   7. 双层融合 (Dual-Layer Blending)：底部低DOPU区域使用稍大核的增强层与基础层融合，保留颗粒感同时抑制杂乱方向
%   8. 统一输出自适应滤波逻辑，缩小核上限以避免底部过度平滑（视觉上更自然，与上方纹理一致）

    %% 参数默认值与初始化
    if nargin < 19 || isempty(bottomDopuThreshold)
        bottomDopuThreshold = 0.5; % 默认底层DOPU阈值
    end
    if nargin < 18 || isempty(adaptiveFilterBottomDepth)
        adaptiveFilterBottomDepth = 120;  % 默认DOPU自适应滤波仅对底部120层生效（0表示全局生效）
    end
    if nargin < 17 || isempty(bottomLayerDepth)
        bottomLayerDepth = 80;  % 默认底层深度为80层
    end
    if nargin < 16 || isempty(enableBottomLayerPhaseReduction)
        enableBottomLayerPhaseReduction = false;  % 默认禁用底层相位延迟减小
    end
    if nargin < 15 || isempty(outputDopuThreshold)
        outputDopuThreshold = 0.4;  % 默认输出滤波DOPU阈值
    end
    if nargin < 14 || isempty(kRU_output)
        kRU_output = 21;  % 默认上限
    end
    if nargin < 13 || isempty(kRL_output)
        kRL_output = 13;   % 默认下限（通常设为h2核大小）
    end
    if nargin < 12 || isempty(enableOutputAdaptive)
        enableOutputAdaptive = false;  % 默认禁用自适应滤波
    end
    if nargin < 11 || isempty(dopuThreshold)
        dopuThreshold = 0.4;  % 默认阈值（用户可传入自定义值）
    end
    if nargin < 10, enableDopuPhaseSupp = true; end
    if nargin < 9 || isempty(dopuMap), dopuMap = []; end
    
    [nZ, nX] = size(Qm);
    
    % 常量定义
    stLayer = 1;
    DTheta = 1.0675 * 2 * pi;
    dTheta = DTheta / 320;
    LABG = [1 0 0];
    
    % 修改1: 移除原有的 median(validDopu) * 1.2 自动阈值计算
    % 改为优先使用函数传入的 dopuThreshold 参数，否则默认 0.4
    dopuLowThresh = dopuThreshold;
    dopuMinScale = 0.1;

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
            
            % 修改2: 直接信任并读取预处理好的 DopuF(i, j)，移除循环内冗余计算
            % 假设外部传入的 DOPU Map 是准确的，不再重复计算
            winDOPU = DopuF(i, j);
            
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
            
            % 修改3: 彻底重写光轴计算分支，基于 DOPU 的自适应分级拟合
            % 
            % 分支 A (winDOPU > dopuLowThresh): 高信噪比区域
            %   - 使用标准多点拟合 (multiPointFitting)
            %   - 拟合点数 = 5，以保留最高的轴向分辨率
            %   - 不应用切面圆半径约束
            %
            % 分支 B (winDOPU <= dopuLowThresh): 噪声/深层区域
            %   - 使用最大窗口拟合 (nFit = Avnum)，利用更多点来压制噪声
            %   - 启用 planeFitConstrained，r_min = 0.6
            %   - 防止光轴随机跳变
            
            if winDOPU > dopuLowThresh
                % === 分支 A: 高 DOPU 区域（高信噪比） ===
                nFit = min(3, Avnum);  % 使用 3 个点拟合（小窗口，高分辨率）
                
                % 取窗口头部的 nFit 个点
                S_use_mod = Stokes_DeMod(1:nFit, :);
                S_use_dc  = Stokes_DeDC(1:nFit, :);
                
                % 使用标准多点拟合，不应用半径约束
                [B_raw_vec, ~] = multiPointFitting(S_use_mod);
                [B_dc_vec, ~]  = multiPointFitting(S_use_dc);
            else
                % === 分支 B: 低 DOPU 区域（低信噪比） ===
                nFit = Avnum;  % 使用所有点拟合（大窗口，压制噪声）
                
                % 取窗口所有点
                S_use_mod = Stokes_DeMod;
                S_use_dc  = Stokes_DeDC;
                
                % 使用切面圆约束拟合，r_min = 0.6（强制半径法）
                r_min = 0.6;
                [B_raw_vec, ~, ~] = planeFitConstrained(S_use_mod, r_min);
                [B_dc_vec, ~, ~]  = planeFitConstrained(S_use_dc, r_min);
            end

            % 使用约束后的向量计算相位和原始累积 LA
            phR_raw(i, j) = calcPhase(S_use_mod, B_raw_vec);
            cumLA_raw(i, j, :) = -B_raw_vec;

            % 基本的窗口跨度度量（使用拟合的子集）
            P1 = S_use_mod(1, :);
            P3 = S_use_mod(end, :);
            T3_nor = norm(P3 - P1);

            % B_dedc 用于后续 DeltaB 逻辑
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

    %% Step 3.6: 底层相位延迟减小 (修正为绝对深度)
    % 在底层深度区域，对低DOPU像素降低相位延迟50%，并进行平滑处理
    if enableBottomLayerPhaseReduction
        % 使用 nZ (原始图像高度) 计算绝对深度阈值
        absDepthThresh = max(1, nZ - bottomLayerDepth + 1);
        
        % 构建基于绝对深度的掩码
        % 展平空间中：绝对深度 = 展平索引 + 表面位置 - 1
        % 条件: FlattenZ + Top - 1 >= absDepthThresh
        % 即: FlattenZ >= absDepthThresh - Top + 1
        maskBottomRegion = false(output_depth, nX);
        for j = 1:nX
            top = test_seg_top(j);
            startIdx = max(1, absDepthThresh - top + 1);
            if startIdx <= output_depth
                maskBottomRegion(startIdx:end, j) = true;
            end
        end
        
        % 结合DOPU阈值（使用专门的底层DOPU阈值，来自配置）
        maskBottomLowDopu = maskBottomRegion & (DopuF(1:output_depth, :) > 0) & (DopuF(1:output_depth, :) < bottomDopuThreshold);
        
        % 定义目标衰减强度 (0.1~0.2 之间效果最好，0.15 为青色，0.1 为蓝色)
        targetScale = 0.15;  % 可根据需要微调（例如 0.15 推荐）
        % 创建减小因子图（底层低DOPU区域使用 targetScale 衰减）
        reductionFactor = ones(output_depth, nX);
        reductionFactor(maskBottomLowDopu) = targetScale;
        
        % 对减小因子进行高斯平滑，避免突兀的边界
        smoothKernel = fspecial('gaussian', [7 7], 2.0);
        reductionFactor_smooth = imfilter(reductionFactor, smoothKernel, 'replicate');
        
        % 将平滑结果限幅到 [targetScale, 1]，以避免平滑引入超出预期的极端值
        reductionFactor_smooth = min(max(reductionFactor_smooth, targetScale), 1);
        
        % 应用平滑后的减小因子到相位延迟（保留底层纹理但降低亮度）
        phR = phR .* reductionFactor_smooth;
        phR_raw = phR_raw .* reductionFactor_smooth;
        phR_rmBG = phR_rmBG .* reductionFactor_smooth;
    end

    %% Step 4: 空间滤波 (双层融合策略)
    % 基础层 (Layer 1): 全图使用标准小核 h2 滤波，保持整体纹理一致性
    % 增强层 (Layer 2): 使用稍大的核 h_strong 滤波，仅在底部低DOPU区域柔和融合
    % 融合公式: final = base .* (1 - mask) + strong .* mask
    
    % --- Layer 1: 基础层滤波 (全图 h2) ---
    [cumLA_bg_base] = filterVectorField(cumLA_bg, h2);
    [cumLA_raw_base] = filterVectorField(cumLA_raw, h2);
    phR_base      = imfilter(phR, h2, 'replicate');
    phR_rmBG_base = imfilter(phR_rmBG, h2, 'replicate');
    phR_raw_base  = imfilter(phR_raw, h2, 'replicate');
    
    if enableOutputAdaptive
        % --- Layer 2: 增强层滤波 (全图 h_strong) ---
        % h_strong 尺寸为 h2 的 1.5~2.0 倍，不要太大以避免塑料感
        h2_size = size(h2, 1);
        h_strong_size = round(h2_size * 1.8);  % 约 1.8 倍
        if mod(h_strong_size, 2) == 0
            h_strong_size = h_strong_size + 1;  % 确保为奇数
        end
        h_strong_sigma = h_strong_size / 5.5;  % sigma 适度增加
        h_strong = fspecial('gaussian', [h_strong_size h_strong_size], h_strong_sigma);
        
        [cumLA_bg_strong] = filterVectorField(cumLA_bg, h_strong);
        [cumLA_raw_strong] = filterVectorField(cumLA_raw, h_strong);
        phR_strong      = imfilter(phR, h_strong, 'replicate');
        phR_rmBG_strong = imfilter(phR_rmBG, h_strong, 'replicate');
        phR_raw_strong  = imfilter(phR_raw, h_strong, 'replicate');
        
        % --- 融合掩码 (Blending Mask) ---
        % 复用 Step 3.6 的 maskBottomLowDopu 逻辑（底部 + 低DOPU）
        % 如果 enableBottomLayerPhaseReduction 未启用，则在此处重新构建掩码
        DopuF_cropped = DopuF(1:output_depth, :);
        
        if ~exist('maskBottomLowDopu', 'var') || isempty(maskBottomLowDopu)
            % 重新构建底部区域掩码
            absDepthThresh_blend = max(1, nZ - adaptiveFilterBottomDepth + 1);
            maskBottomRegion_blend = false(output_depth, nX);
            for j = 1:nX
                top = test_seg_top(j);
                startIdx = max(1, absDepthThresh_blend - top + 1);
                if startIdx <= output_depth
                    maskBottomRegion_blend(startIdx:end, j) = true;
                end
            end
            % 结合DOPU阈值
            maskBottomLowDopu_blend = maskBottomRegion_blend & (DopuF_cropped > 0) & (DopuF_cropped < outputDopuThreshold);
        else
            maskBottomLowDopu_blend = maskBottomLowDopu;
        end
        
        % --- 掩码羽化 (Feathering) ---
        % 使用高斯模糊对掩码进行软化，避免融合边界生硬
        featherKernel = fspecial('gaussian', [15 15], 4.0);
        blendMask = double(maskBottomLowDopu_blend);
        blendMask = imfilter(blendMask, featherKernel, 'replicate');
        blendMask = min(max(blendMask, 0), 1);  % 确保在 [0, 1] 范围内
        
        % --- 双层融合 ---
        % 相位延迟的融合（标量）
        phR_gF      = phR_base .* (1 - blendMask) + phR_strong .* blendMask;
        phR_rmBG_gF = phR_rmBG_base .* (1 - blendMask) + phR_rmBG_strong .* blendMask;
        phR_raw_gF  = phR_raw_base .* (1 - blendMask) + phR_raw_strong .* blendMask;
        
        % 光轴向量的融合（需要重新归一化）
        blendMask3 = repmat(blendMask, [1 1 3]);  % 扩展到3通道
        cumLA_bg_blended = cumLA_bg_base .* (1 - blendMask3) + cumLA_bg_strong .* blendMask3;
        cumLA_raw_blended = cumLA_raw_base .* (1 - blendMask3) + cumLA_raw_strong .* blendMask3;
        
        % 重新归一化光轴向量
        cumLA_bg_gF = normalizeVectorField(cumLA_bg_blended);
        cumLA_raw_gF = normalizeVectorField(cumLA_raw_blended);
    else
        % 非自适应模式：直接使用基础层结果
        cumLA_bg_gF = cumLA_bg_base;
        cumLA_raw_gF = cumLA_raw_base;
        phR_gF      = phR_base;
        phR_rmBG_gF = phR_rmBG_base;
        phR_raw_gF  = phR_raw_base;
    end

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

%% 向量场归一化（用于融合后重新归一化）
function [V_out] = normalizeVectorField(V_in)
    V1 = V_in(:,:,1);
    V2 = V_in(:,:,2);
    V3 = V_in(:,:,3);
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