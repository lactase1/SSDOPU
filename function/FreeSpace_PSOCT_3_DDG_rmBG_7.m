% ============================================================================
% FreeSpace_PSOCT_3_DDG_rmBG_7.m - 基于DDG的双折射计算函数
% ============================================================================
%
% 功能描述:
%   使用差分梯度分解(DDG)方法从PS-OCT数据计算双折射参数
%   包括光轴方向(Optical Axis)和相位延迟(Phase Retardation)
%   支持背景双折射校正和深度分辨率计算
%
% 主要特性:
%   - 高斯滤波代替传统中值和平滑滤波
%   - QUV小尺度滤波(h1)和光轴大尺度滤波(h2)
%   - 每个A-scan的DC分量去除和向量标准化
%   - 使用滤波后光轴计算旋转矩阵
%   - 输出累积光轴 - 2024/08/25
%   - 基于信噪比(DOPU)动态调整拟合圆盘直径，避免噪声导致的相位延迟高估
%
% 算法特性:
%   - 使用DDG方法计算光轴和相位延迟
%   - 以T3(P3-P1的范数)作为背景和偏振信号的区分标记
%   - 使用norm(P3-P1)校正背景效应
%   - 输出去除背景后的相位延迟
%   - 计算深度分辨率光轴需要包含背景诱导相位的相位延迟 - 2024/08/31
%   - 基于信噪比(DOPU)动态调整拟合圆盘直径，避免噪声导致的相位延迟高估
%
% 历史更新:
%   - 2024/09/02: 测试视网膜成像
%   - 使用超过3个数据点计算光轴并平均相位延迟
%   - 支持体积数据处理
%
% 输入参数:
%   Qm, Um, Vm: Q、U、V Stokes参数矩阵
%   Stru: 结构矩阵
%   test_seg_top: 分割顶部位置
%   h1, h2: 滤波核
%   Avnum: 最大平均层数(窗口大小)，决定在每个深度位置最多使用多少个连续层来计算双折射
%   dopu_splitSpec_M: 分光谱DOPU矩阵，用于评估信噪比并动态调整拟合圆盘直径
%
% 输出参数:
%   LA, PhR, cumLA: 背景校正后的深度分辨率光轴、相位延迟、累积光轴
%   LA_raw, PhR_raw, cumLA_raw: 原始结果(无背景校正)
% ============================================================================

function [LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = FreeSpace_PSOCT_3_DDG_rmBG_7(Qm,Um,Vm,Stru,test_seg_top,h1,h2,Avnum)
showimg = 0;


%% 数据展平处理 (Data Flattening)
% 将基于曲面分割的数据转换为平面格式
% OCT数据通常沿着曲面采集，从组织表面开始而不是图像顶部
% 将每个B-scan从分割点开始的数据"展平"到统一的高度起始点
% 这样所有B-scan都从"相对深度=1"开始，便于后续统一处理

% 计算展平后的最大长度：确保数组足够大
max_flattened_length = max(size(Stru,1) - min(test_seg_top) + 1, size(Stru,1));

% 预分配展平数组
StruF = zeros(max_flattened_length, size(Stru,2));
QmF = zeros(max_flattened_length, size(Qm,2));
UmF = zeros(max_flattened_length, size(Um,2));
VmF = zeros(max_flattened_length, size(Vm,2));

% 实际填充数据
for j=1:size(Stru,2)
   start_idx = max(1, test_seg_top(j));  % 确保起始索引至少为1
   flattened_length = size(Stru,1) - start_idx + 1;
   StruF(1:flattened_length, j) = Stru(start_idx:size(Stru,1), j);
   QmF(1:flattened_length, j) = Qm(start_idx:size(Qm,1), j);
   UmF(1:flattened_length, j) = Um(start_idx:size(Um,1), j);
   VmF(1:flattened_length, j) = Vm(start_idx:size(Vm,1), j);
end


%% QUV标准化和预处理 (QUV Normalization and Preprocessing)
% 对QUV数据进行标准化处理，便于后续向量计算
% 将QUV向量归一化到单位球面上
QF=QmF;UF=UmF;VF=VmF;

QF1=QF./sqrt(QF.^2+UF.^2+VF.^2); % 向量标准化
UF1=UF./sqrt(QF.^2+UF.^2+VF.^2);
VF1=VF./sqrt(QF.^2+UF.^2+VF.^2);

QF1(isnan(QF1)) = 0;    UF1(isnan(UF1)) = 0;    VF1(isnan(VF1)) = 0;

SmapF=cat(3,QF1,UF1,VF1); % 合并为3D矩阵 [深度×X×3]
% figure(11),imagesc(SmapF);

SmapF1=smoothdata(SmapF, 1, 'movmean', 15); % 沿深度方向进行移动平均平滑，减少噪声


%% DDG双折射计算 (DDG Birefringence Calculation)
% 使用Differential Decomposition of the Gradient方法计算双折射
% 从多层QUV数据点拟合双折射轴和相位延迟
stLayer =1;
DTheta = 1.0675*2*pi; % 背景最大相位延迟估计 [rad] (50cm路径导致)
dTheta = DTheta/320; % 每像素相位延迟 [rad/px]
LABG = [1 0 0]; % 背景双折射轴方向(沿X轴)
% Avnum = 3;
phd=zeros(size(SmapF1,1),size(SmapF1,2));

% 对每个A-scan进行计算
for j=1:size(SmapF1,2)  % X方向循环
    dirCheck = 1;
    for i=1:size(SmapF1,1)-3  % 深度方向循环(最小Avnum=3，避免越界)
        % 使用固定Avnum
        dynamic_Avnum = Avnum;
        
        % 确保不会越界
        max_possible_i = size(SmapF1,1) - dynamic_Avnum;
        if i > max_possible_i
            dynamic_Avnum = max(3, size(SmapF1,1) - i);
        end
        
        % 迭代计算：如果轨迹半径太小，增加Avnum重新计算
        min_radius_threshold = 0.05;  % 最小轨迹半径阈值
        max_iterations = 3;  % 最大迭代次数
        iteration = 1;
        valid_calculation = false;
        
        while ~valid_calculation && iteration <= max_iterations
            planeData=squeeze(SmapF1(i:i+dynamic_Avnum-1,j,:)); % 提取dynamic_Avnum个连续层的QUV数据
            x1=squeeze(planeData(:,1)); y1=squeeze(planeData(:,2)); z1=squeeze(planeData(:,3));
            nx=sum(sum(x1==0));
            if x1==0|nx>0 % 数据有效性检查
                cumLA_bg(i,j,1:3)=0;
                cumLA_raw(i,j,1:3)=0;
                phR(i,j)=0;
                phR_rmBG(i,j) = 0;
                phR_raw(i,j) = 0;
                phd(i,j)=0;
                rotationVector(i,j,1:3)=0;
                valid_calculation = true;  % 无数据，标记为有效
                break;
            end
            
            % 背景校正计算
            ps_bg_rdu = [];ps_bg_rm = [];
            for jL =1:dynamic_Avnum
                inpp = planeData(jL,:);
                outp = rotationVectorToMatrix(-1*(i+test_seg_top(j)+stLayer-0)*dTheta*LABG)*inpp';
                outp2 = rotationVectorToMatrix(-1*(i+test_seg_top(j)+stLayer+jL-2)*dTheta*LABG)*inpp';
                ps_bg_rdu = [ps_bg_rdu;outp'];
                ps_bg_rm = [ps_bg_rm;outp2'];
            end
            
            % 选择三个代表性点
            P1  = ps_bg_rm(1,:); 
            P2 = ps_bg_rm(max(1, round((1+dynamic_Avnum)/2)),:); 
            P3 = ps_bg_rm(dynamic_Avnum,:);
            P1p = ps_bg_rdu(1,:);
            P2p = ps_bg_rdu(min(size(ps_bg_rdu,1), 2),:);
            P3p = ps_bg_rdu(min(size(ps_bg_rdu,1), 3),:);
            T1p = P2p - P1p; T2p = P3p-P2p; T1 = P2-P1;T2 = P3-P2;
            Tb1 = P1 -P1p; Tb2 = P3-P3p;
            
            % 计算轨迹半径（P3-P1的范数）
            trajectory_radius = norm(P3-P1);
            
            % 如果轨迹半径太小，增加Avnum重新计算
            if trajectory_radius < min_radius_threshold && dynamic_Avnum < Avnum && iteration < max_iterations
                dynamic_Avnum = min(dynamic_Avnum + 2, Avnum);  % 增加Avnum
                % 重新检查越界
                max_possible_i_new = size(SmapF1,1) - dynamic_Avnum;
                if i > max_possible_i_new
                    dynamic_Avnum = max(3, size(SmapF1,1) - i);
                end
                iteration = iteration + 1;
                continue;  % 重新计算
            end
            
            % 轨迹半径足够大，继续正常计算
            valid_calculation = true;
            
            % DDG计算: 准备点和向量用于双折射计算
            % 选择三个代表性点: 起点、中点、终点
            P1  = ps_bg_rm(1,:); 
            P2 = ps_bg_rm(max(1, round((1+dynamic_Avnum)/2)),:); 
            P3 = ps_bg_rm(dynamic_Avnum,:);
            P1p = ps_bg_rdu(1,:);
            P2p = ps_bg_rdu(min(size(ps_bg_rdu,1), 2),:);
            P3p = ps_bg_rdu(min(size(ps_bg_rdu,1), 3),:);
            T1p = P2p - P1p; T2p = P3p-P2p; T1 = P2-P1;T2 = P3-P2;
            Tb1 = P1 -P1p; Tb2 = P3-P3p;
            
            % 从背景校正点直接计算双折射和相位延迟 [bgrm DDG]
            % 双折射轴计算: 两个切向量的叉积给出垂直于轨迹平面的法向量
            B_bg_rm = cross(T1,T2); B_bg_rm = B_bg_rm/sqrt(sum(B_bg_rm.^2));
            % 使用法向量计算相位延迟
            [Subphr_bg_rm_Ts] = ddgTsPhr(ps_bg_rm,B_bg_rm);
            phR_raw(i,j) = acos(median(Subphr_bg_rm_Ts));% avg acos first, then acos
            cumLA_raw(i,j,1:3) = -B_bg_rm;
            
            % 背景校正DDG: 先计算Bp，再用deltaB校正
            % 原理: 区分背景双折射和组织双折射
            Bp = cross(T1p,T2p); 
            T3_nor = norm(P3-P1); % P3-P1的范数，用于区分背景和信号
            % 根据T3_nor判断是否需要背景校正
            if T3_nor < 0.065%0.7/0.12  % 阈值判断
                deltaB = [0 0 0]; % 无需校正
                T1 = T1p; T2 = T2p;
                if dirCheck && (Bp(1) > 0) %% check LA direction
                    Bp = -Bp;
                end
            else
                dirCheck = 0;
                deltaB = cross(T1p,Tb2) + cross(T2p,Tb1) - cross(Tb1,Tb2);% 背景校正计算
            end
            B = Bp + deltaB; % 校正后的双折射轴
            B = B/sqrt(sum(B.^2)); % 标准化
            
            T3_nors(i,j) = T3_nor;
            if 1  % 总是执行的分支
                if T3_nor < 0.065
                    phR(i,j) = acos(median(ddgTsPhr(ps_bg_rdu,B)));
                    phR_rmBG(i,j) = max([phR(i,j) - 1*dTheta, 0]); % 去除背景相位
                    phR_rmBG(i,j) = min([phR(i,j), 0.2]); % 上限限制
                else
                    phR(i,j) = acos(median(ddgTsPhr(ps_bg_rm,B)));
                    phR_rmBG(i,j) = phR(i,j);
                end
                cumLA_bg(i,j,1:3)=-B; % 背景校正后的累积双折射
            else  % 从不执行的分支
                N1 = cross(B,T1); N2 = cross(B,T2);% normal vector
                phr = acos(dot(N1,N2)/sqrt(sum(N1.^2)*sum(N2.^2)));% phase retardation
                phR(i,j)=phr;
                phr = max([phr-1*dTheta, 0]);
                cumLA_bg(i,j,1:3)=-B;
                phR_rmBG(i,j) = phr; %% rm BG 
                Bps(i,j,1:3) = Bp/sqrt(sum(Bp.^2));
                T3_nors(i,j) = T3_nor;
            end
            
        end
    end
end


%% 后处理: 对计算结果进行滤波和深度分辨率转换
% 光轴滤波: 减少噪声，提高双折射轴估计的稳定性
[cumLA_bg_gF] = LAgFfilt(cumLA_bg,h2);
[cumLA_raw_gF] = LAgFfilt(cumLA_raw,h2);
if showimg
figure;imshow(cumLA_bg_gF*0.5+0.5,[])
figure;imshow(cumLA_raw_gF*0.5+0.5,[])
figure;imshow(phR_rmBG,[0 0.25])
figure;imshow(phR_raw,[0 0.25])
end

%% 相位滤波: 减少相位延迟估计的噪声
    phR_gF=imfilter(phR,h2,'replicate');%% with BG
    phR_rmBG_gF=imfilter(phR_rmBG,h2,'replicate');%% with BG
    phR_raw_gF = imfilter(phR_raw,h2,'replicate');%% raw w/o BG
    phrRg = [0 0.3];
    if showimg
    figure;imshowpair(mat2gray(phR_rmBG,phrRg),mat2gray(phR_rmBG_gF,phrRg),'montage')
     figure;imshowpair(mat2gray(phR_raw,phrRg),mat2gray(phR_raw_gF,phrRg),'montage')
    end

%% 计算旋转矩阵: 为深度分辨率转换计算初始旋转矩阵
% for i=1:size(SmapF1,1)-Avnum
    for j=1:size(SmapF1,2)
        ax_rot(1,:)=cumLA_bg_gF(1,j,:);
        rotationVector(1,j,:) = -phR_gF(1,j)/2 * ax_rot;
        ax_rot(1,:)=cumLA_raw_gF(1,j,:);
        rotationVector_raw(1,j,:) = -phR_raw_gF(1,j)/2 * ax_rot;
    end
% end
% loaxis22=loaxis3.*(phR);
% loaxis22=cumLA_bg_gF;


%% 深度分辨率转换: 将累积光轴转换为深度分辨率光轴
% 从累积双折射得到每层深度处的局部双折射轴
% 通过旋转矩阵逆运算，将累积效应转换为局部效应
[bfloaxis3D] = drLA(cumLA_bg_gF,-phR_gF,rotationVector);
[bfloaxis3D_raw] = drLA(cumLA_raw_gF,-phR_raw_gF,rotationVector_raw);
% [bfloaxis3D] = drLA(cumLA_bg_gF.*phR_gF,-phR_gF,rotationVector);
% [bfloaxis3D] = drLA(cumLA_bg_gF.*phR_rmBG_gF,-phR_rmBG_gF,rotationVector);
% [bfloaxis3D_raw] = drLA(cumLA_raw_gF.*phR_raw_gF,-phR_raw_gF,rotationVector_raw);
if showimg
    figure;imshowpair(cumLA_bg_gF*0.5+0.5,bfloaxis3D*0.5+0.5,'montage')
    figure;imshowpair(cumLA_raw_gF*0.5+0.5,bfloaxis3D_raw*0.5+0.5,'montage')
    phrRg = [0 0.3];
    figure;imshowpair(mat2gray(phR_rmBG_gF,phrRg),mat2gray(phR_raw_gF,phrRg),'montage')
end

%%
% with BG
bfphrr=phR_rmBG_gF;
bfloaxis3D(isnan(bfloaxis3D)) = 0;
bfphrr(isnan(bfphrr)) = 0;
% w/o BG
bfphrr_raw=phR_raw_gF;
bfloaxis3D_raw(isnan(bfloaxis3D_raw)) = 0;
bfphrr_raw(isnan(bfphrr_raw)) = 0;

% 计算正确的输出数组大小：基于原始输入数据大小减去Avnum
expected_output_size = size(Qm, 1) - Avnum;

% 初始化临时数组，保持原始大小用于放置逻辑
temp_size = size(bfloaxis3D, 1);
temp_LA = zeros(temp_size, size(Qm, 2), 3);
temp_PhR = zeros(temp_size, size(Qm, 2));
temp_cumLA = zeros(temp_size, size(Qm, 2), 3);
temp_LA_raw = zeros(temp_size, size(Qm, 2), 3);
temp_PhR_raw = zeros(temp_size, size(Qm, 2));
temp_cumLA_raw = zeros(temp_size, size(Qm, 2), 3);

% restore curved results - 使用原始放置逻辑
for j=1:size(Qm,2)
    toInd = test_seg_top(j):min(temp_size, test_seg_top(j) + size(bfloaxis3D, 1) - 1);
    fromInd = 1:length(toInd);
    if ~isempty(toInd) && ~isempty(fromInd)
        temp_LA(toInd, j, :) = bfloaxis3D(fromInd, j, :);
        temp_PhR(toInd, j) = bfphrr(fromInd, j);
        temp_cumLA(toInd, j, :) = cumLA_bg_gF(fromInd, j, :);

        temp_LA_raw(toInd, j, :) = bfloaxis3D_raw(fromInd, j, :);
        temp_PhR_raw(toInd, j) = bfphrr_raw(fromInd, j);
        temp_cumLA_raw(toInd, j, :) = cumLA_raw_gF(fromInd, j, :);
    end
end

% 截断到正确的大小 - 取前 expected_output_size 行
LA = temp_LA(1:expected_output_size, :, :);
PhR = temp_PhR(1:expected_output_size, :);
cumLA = temp_cumLA(1:expected_output_size, :, :);
LA_raw = temp_LA_raw(1:expected_output_size, :, :);
PhR_raw = temp_PhR_raw(1:expected_output_size, :);
cumLA_raw = temp_cumLA_raw(1:expected_output_size, :, :);
end

%% 使用切向量计算相位延迟: 计算相邻点之间的相位延迟
% 输入: ps_bg_rm - 背景校正后的3D点序列, B_bg_rm - 双折射轴
% 输出: Subphr_bg_rm_Ts - 相位延迟余弦值数组
function [Subphr_bg_rm_Ts] = ddgTsPhr(ps_bg_rm,B_bg_rm)
% 3D points
% LA
    Avnum = size(ps_bg_rm,1);
    for k = 1:Avnum-2
        T1_bg_rm = ps_bg_rm(k+1,:)-ps_bg_rm(k,:);
        T2_bg_rm = ps_bg_rm(k+2,:)-ps_bg_rm(k+1,:);
        N1_bg_rm = cross(B_bg_rm,T1_bg_rm);N2_bg_rm = cross(B_bg_rm,T2_bg_rm);% normal vector
        Subphr_bg_rm_Ts(k) = dot(N1_bg_rm,N2_bg_rm)/(sqrt(sum(N1_bg_rm.^2)*sum(N2_bg_rm.^2)));
    end
    Subphr_bg_rm_Ts(Subphr_bg_rm_Ts > 1) = 1; Subphr_bg_rm_Ts(Subphr_bg_rm_Ts < -1) = -1;
end

%% 使用点积计算相位延迟: 计算相邻点在垂直于光轴平面上的投影
function [Subphr_bg_rm_Ps] = ddgPsPhr(ps_bg_rm,B_bg_rm)
% 3D points
% LA
    Avnum = size(ps_bg_rm,1);
    for k = 1:Avnum-1
        P1t_bg_rm = cross(ps_bg_rm(k,:),B_bg_rm);
        P2t_bg_rm = cross(ps_bg_rm(k+1,:),B_bg_rm);
        Subphr_bg_rm_Ps(k) = dot(P1t_bg_rm,P2t_bg_rm)/(sqrt(sum(P1t_bg_rm.^2)*sum(P2t_bg_rm.^2)));
    end
    Subphr_bg_rm_Ps(Subphr_bg_rm_Ps > 1) = 1; Subphr_bg_rm_Ps(Subphr_bg_rm_Ps < -1) = -1;
end

%% 光轴滤波函数: 对3D光轴向量进行高斯滤波并重新标准化
function [loaxis3] = LAgFfilt(loaxis2,h2)
    Temp_ax1(:,:)=loaxis2(:,:,1);
    Temp_ax2(:,:)=loaxis2(:,:,2);
    Temp_ax3(:,:)=loaxis2(:,:,3);

    Temp_ax1=imfilter(Temp_ax1,h2,'replicate');
    Temp_ax2=imfilter(Temp_ax2,h2,'replicate');
    Temp_ax3=imfilter(Temp_ax3,h2,'replicate');

    Temp_ax1=Temp_ax1./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
    Temp_ax2=Temp_ax2./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
    Temp_ax3=Temp_ax3./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);

    Temp_ax1(isnan(Temp_ax1)) = 0;    Temp_ax2(isnan(Temp_ax2)) = 0;    Temp_ax3(isnan(Temp_ax3)) = 0;
    loaxis3=cat(3,Temp_ax1,Temp_ax2,Temp_ax3);
end

%% 深度分辨率转换函数: 将累积光轴转换为深度分辨率光轴
% 从组织表面开始，逐步"展开"累积的双折射效应
function [axiss2] = drLA(loaxis22,phR_gF,rotationVector2)

for j=1:size(loaxis22,2)
    for i=1:size(loaxis22,1)
        vector1=squeeze(rotationVector2(1,j,:));
        a=[loaxis22(i,j,1),loaxis22(i,j,2),loaxis22(i,j,3)];
        if i==1
            rotationMatrix1=rotationVectorToMatrix(vector1);
            axis2(i,j,1:3)=a;
        else
            axis2(i,j,1:3)=rotationMatrix1*a';  %当前光轴a绕着上一个光轴ax_rot旋转-phR(i,j)/2角度。
            Qa=axis2(i,j,1);    Ua=axis2(i,j,2);   Va=axis2(i,j,3);
            d1= -phR_gF(i,j)/2.0*[Qa,Ua,Va]; %[Qa,Ua,Va]是旋转完的光轴，
            rotationMatrix1 = rotationVectorToMatrix(d1)*rotationMatrix1;
%                        rotationMatrix1 = rotationVectorToMatrix(d1);
            % post-process to visualization purpose
            Q2=Qa/sqrt(Qa^2+Ua^2+Va^2);  U2=Ua/sqrt(Qa^2+Ua^2+Va^2);  V2=Va/sqrt(Qa^2+Ua^2+Va^2); 
            V2(isnan(V2)) = 0;
%              Q2=Qa/sqrt(Qa^2+Ua^2);  U2=Ua/sqrt(Qa^2+Ua^2);   
%              Q2(isnan(Q2)) = 0;U2(isnan(U2)) = 0;
            axiss2(i,j,:)=cat(3,Q2,U2,V2);
        end
    end
end

end


