% ============================================================================
% FreeSpace_PSOCT_3_DDG_rmBG_7.m - PS-OCT 极化参数计算函数
%
% 功能描述:
%   该函数实现基于DDG (Double Difference Geometry) 算法的PS-OCT极化参数计算，
%   包括背景去除和校正处理。函数返回两种处理版本的结果：完整背景校正版本
%   和部分背景处理版本，用于对比分析。
%
% 输入参数:
%   Qm, Um, Vm: 极化状态的Stokes参数矩阵 (Q, U, V分量)
%   Stru: 结构图像 (OCT强度图像)
%   test_seg_top: 每个A-line的表面分割位置 (用于坐标变换)
%   h1: 小尺度滤波器 (用于QUV滤波，当前代码中未使用)
%   h2: 大尺度滤波器 (用于光轴和相位滤波)
%   Avnum: 每个深度位置使用的相邻层数 (用于DDG计算)
%
% 输出参数:
%   ============================================================================
%   | 完整背景校正版本 (推荐用于最终分析)                                    |
%   ============================================================================
%   LA:     深度分辨双折射轴 (Linear Axis)
%           - 基于完整背景校正的cumLA_bg_gF计算
%           - 通过drLA函数进行深度旋转，得到每个深度位置的光轴方向
%           - 轴方向经过deltaB校正，相位用于旋转的phR_gF保留BG诱导项
%
%   PhR:    相位延迟 (Phase Retardation)
%           - 经过完整背景校正的相位延迟 (phR_rmBG_gF)
%           - 根据T3_nor条件进行背景校正处理
%
%   cumLA:  累积双折射轴 (Cumulative Linear Axis)
%           - 经过完整背景校正的累积双折射轴 (cumLA_bg_gF)
%           - 通过LAgFfilt进行光轴滤波
%
%   ============================================================================
%   | 部分背景处理版本 (用于算法对比和诊断)                                |
%   ============================================================================
%   LA_raw: 原始深度分辨双折射轴
%           - 基于部分背景处理的cumLA_raw_gF计算
%           - 只去除DC背景分量，无deltaB复杂校正
%           - 直接从去DC点集计算DDG结果
%
%   PhR_raw: 原始相位延迟
%           - 经过部分背景处理的相位延迟 (phR_raw_gF)
%           - 直接基于去DC点集计算，无完整背景校正
%
%   cumLA_raw: 原始累积双折射轴
%           - 经过部分背景处理的累积双折射轴 (cumLA_raw_gF)
%           - 直接从DDG计算，无deltaB校正
%
% 处理流程概述:
%   1. 数据预处理: 坐标变换和QUV标准化
%   2. SVD光轴拟合: 从QUV数据估计初始光轴
%   3. 背景处理: 去除DC分量并进行复杂背景校正
%   4. DDG计算: 使用双差分几何计算双折射参数
%   5. 后处理: 滤波和深度分辨率恢复
%
% 版本历史:
%   - 基于DDG算法的PS-OCT极化分析
%   - 支持完整背景校正和部分背景处理两种模式
%   - 用于视网膜成像和组织双折射定量分析
%
% ============================================================================

function [LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = FreeSpace_PSOCT_3_DDG_rmBG_7(Qm,Um,Vm,Stru,test_seg_top,h1,h2,Avnum)
showimg = 0;  % 调试显示开关 (0=关闭, 1=开启)

%% 数据预处理: 坐标变换到表面基准坐标系
StruF=Stru*0; QmF=Qm*0; UmF=Um*0; VmF=Vm*0;
for j=1:size(Stru,2)
   StruF(1:size(Stru,1)-test_seg_top(j)+1,j)=Stru(test_seg_top(j):size(Stru,1),j);
   QmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Qm(test_seg_top(j):size(Stru,1),j);
   UmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Um(test_seg_top(j):size(Stru,1),j);
   VmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Vm(test_seg_top(j):size(Stru,1),j);
end

%% QUV数据标准化处理
QF=QmF; UF=UmF; VF=VmF;  % 使用原始数据 (可选: 可在此添加滤波)
QF1=QF./sqrt(QF.^2+UF.^2+VF.^2);  % 归一化Stokes向量
UF1=UF./sqrt(QF.^2+UF.^2+VF.^2);
VF1=VF./sqrt(QF.^2+UF.^2+VF.^2);
QF1(isnan(QF1)) = 0; UF1(isnan(UF1)) = 0; VF1(isnan(VF1)) = 0;

SmapF=cat(3,QF1,UF1,VF1);  % 创建3D极化状态映射
SmapF1=smoothdata(SmapF, 1, 'movmean', 15);  % 深度方向平滑

%% DDG算法参数初始化
stLayer =1;  % 起始层
DTheta = 1.0675*2*pi;  % 背景诱导最大相位延迟估计 (弧度)
dTheta = DTheta/320;   % 每像素相位延迟 (弧度/像素)
LABG = [1 0 0];        % 背景光轴方向 (沿Q轴)

%% 主处理循环: 对每个A-line进行DDG计算
for j=1:size(SmapF1,2)  % 遍历所有A-line
    dirCheck = 1;  % 方向检查标志
    for i=1:size(SmapF1,1)-Avnum  % 遍历深度位置
        planeData=squeeze(SmapF1(i:i+Avnum,j,:));  % 提取Avnum层的数据块
        x1=squeeze(planeData(:,1)); nx=sum(sum(x1==0));

        if x1==0|nx>0  % 跳过无效数据 (Q分量为0)
            cumLA_bg(i,j,1:3)=0; cumLA_raw(i,j,1:3)=0;
            phR(i,j)=0; phR_rmBG(i,j)=0; phR_raw(i,j)=0;
            rotationVector(i,j,1:3)=0;
        else
            %% 步骤1: DC背景去除 - 旋转到去背景坐标系
            ps_bg_rdu = []; ps_bg_rm = [];
            for jL=1:Avnum
                inpp = planeData(jL,:);
                outp = rotationVectorToMatrix(-1*(i+test_seg_top(j)+stLayer-0)*dTheta*LABG)*inpp';
                outp2 = rotationVectorToMatrix(-1*(i+test_seg_top(j)+stLayer+jL-2)*dTheta*LABG)*inpp';
                ps_bg_rdu = [ps_bg_rdu;outp']; ps_bg_rm = [ps_bg_rm;outp2'];
            end

            %% 步骤2: 准备DDG几何元素 (起点/中点/终点)
            P1 = ps_bg_rm(1,:); P2 = ps_bg_rm(round((1+Avnum)/2),:); P3 = ps_bg_rm(Avnum,:);
            P1p = ps_bg_rdu(1,:); P2p = ps_bg_rdu(2,:); P3p = ps_bg_rdu(3,:);
            T1p = P2p - P1p; T2p = P3p-P2p; T1 = P2-P1; T2 = P3-P2;
            Tb1 = P1 -P1p; Tb2 = P3-P3p;

            %% 步骤3: LA_raw计算 (部分背景处理 - 只去DC)
            B_bg_rm = cross(T1,T2); B_bg_rm = B_bg_rm/sqrt(sum(B_bg_rm.^2));
            [Subphr_bg_rm_Ts] = ddgTsPhr(ps_bg_rm,B_bg_rm);
            phR_raw(i,j) = acos(median(Subphr_bg_rm_Ts));
            cumLA_raw(i,j,1:3) = -B_bg_rm;

            %% 步骤4: LA计算 (完整背景校正 - DC去除 + deltaB校正)
            Bp = cross(T1p,T2p); T3_nor = norm(P3-P1);

            % 背景校正逻辑: 根据T3_nor决定校正策略
            if T3_nor < 0.065  % 低噪声情况
                deltaB = [0 0 0]; T1 = T1p; T2 = T2p;
                if dirCheck && (Bp(1) > 0), Bp = -Bp; end  % 方向校正
            else  % 高噪声情况 - 应用复杂deltaB校正
                dirCheck = 0;
                deltaB = cross(T1p,Tb2) + cross(T2p,Tb1) - cross(Tb1,Tb2);
            end

            B = Bp + deltaB; B = B/sqrt(sum(B.^2));  % 应用校正并归一化
            T3_nors(i,j) = T3_nor;

            % 相位计算和背景校正
            if T3_nor < 0.065
                phR(i,j) = acos(median(ddgTsPhr(ps_bg_rdu,B)));
                phR_rmBG(i,j) = max([phR(i,j) - 1*dTheta, 0]);
                phR_rmBG(i,j) = min([phR(i,j), 0.2]);
            else
                phR(i,j) = acos(median(ddgTsPhr(ps_bg_rm,B)));
                phR_rmBG(i,j) = phR(i,j);
            end
            cumLA_bg(i,j,1:3) = -B;
        end
    end
end
%% 后处理: 空间滤波减少噪声
[cumLA_bg_gF] = LAgFfilt(cumLA_bg,h2);   % 完整版本累积轴滤波
[cumLA_raw_gF] = LAgFfilt(cumLA_raw,h2); % 部分版本累积轴滤波

%% 相位滤波处理
phR_gF = imfilter(phR,h2,'replicate');           % 含背景相位 (用于旋转)
phR_rmBG_gF = imfilter(phR_rmBG,h2,'replicate'); % 背景校正相位
phR_raw_gF = imfilter(phR_raw,h2,'replicate');   % 部分处理相位

%% 计算深度旋转矩阵 (为深度分辨率恢复准备)
for j=1:size(SmapF1,2)
    ax_rot(1,:)=cumLA_bg_gF(1,j,:);
    rotationVector(1,j,:) = -phR_gF(1,j)/2 * ax_rot;      % 完整版本旋转向量
    ax_rot(1,:)=cumLA_raw_gF(1,j,:);
    rotationVector_raw(1,j,:) = -phR_raw_gF(1,j)/2 * ax_rot; % 部分版本旋转向量
end

%% 深度分辨率恢复: 通过旋转累积恢复各深度位置的双折射轴
[bfloaxis3D] = drLA(cumLA_bg_gF,-phR_gF,rotationVector);         % 完整背景校正版本
[bfloaxis3D_raw] = drLA(cumLA_raw_gF,-phR_raw_gF,rotationVector_raw); % 部分处理版本

%% 结果格式化: 从平坦化坐标系映射回原始弯曲坐标系
bfphrr = phR_rmBG_gF;      % 完整版本相位
bfphrr_raw = phR_raw_gF;   % 部分版本相位

% 初始化输出数组
LA=bfloaxis3D*0; PhR=bfphrr*0; cumLA=bfloaxis3D*0;
LA_raw=bfloaxis3D*0; PhR_raw=bfphrr*0; cumLA_raw=bfloaxis3D*0;

% 坐标变换: 平坦化坐标系 -> 原始弯曲坐标系
for j=1:size(Stru,2)
    toInd = test_seg_top(j):size(LA,1);           % 原始坐标系目标索引
    fromInd = 1:size(LA,1)-test_seg_top(j)+1;     % 平坦化坐标系源索引

    % 完整背景校正版本映射
    LA(toInd,j,:)=bfloaxis3D(fromInd,j,:);
    PhR(toInd,j)=bfphrr(fromInd,j);
    cumLA(toInd,j,:)=cumLA_bg_gF(fromInd,j,:);

    % 部分背景处理版本映射
    LA_raw(toInd,j,:)=bfloaxis3D_raw(fromInd,j,:);
    PhR_raw(toInd,j)=bfphrr_raw(fromInd,j);
    cumLA_raw(toInd,j,:)=cumLA_raw_gF(fromInd,j,:);
end

%% 输出变量说明:
%% LA: 完整背景校正的双折射轴 (深度×宽度×3)
%% PhR: 完整背景校正的相位延迟 (深度×宽度)
%% cumLA: 完整背景校正的累积双折射轴 (深度×宽度×3)
%% LA_raw: 部分背景校正的双折射轴 (深度×宽度×3)
%% PhR_raw: 部分背景校正的相位延迟 (深度×宽度)
%% cumLA_raw: 部分背景校正的累积双折射轴 (深度×宽度×3)
end

%% ============================================================================
%% 子函数: DDG相位计算 (使用切向量)
%% ============================================================================
%% 函数功能:
%%   使用DDG算法计算相邻点之间的相位延迟
%%   基于双折射轴与切向量的夹角计算相位差
%%
%% 输入参数:
%%   ps_bg_rm: 去背景后的3D点集 (N×3矩阵)
%%   B_bg_rm: 双折射轴向量 (1×3向量)
%%
%% 输出参数:
%%   Subphr_bg_rm_Ts: 相位延迟数组 (弧度)
%%
%% 计算原理:
%%   对于连续的三点P1,P2,P3，计算向量T1=P2-P1, T2=P3-P2
%%   相位延迟 = arccos(dot(N1,N2))，其中N1=cross(B,T1), N2=cross(B,T2)
%% ============================================================================
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

%% ============================================================================
%% 子函数: DDG相位计算 (使用法向量)
%% ============================================================================
%% 函数功能:
%%   使用DDG算法的另一种实现计算相位延迟
%%   基于双折射轴与法向量的夹角计算相位差
%%
%% 输入参数:
%%   ps_bg_rm: 去背景后的3D点集 (N×3矩阵)
%%   B_bg_rm: 双折射轴向量 (1×3向量)
%%
%% 输出参数:
%%   Subphr_bg_rm_Ps: 相位延迟数组 (弧度)
%%
%% 计算原理:
%%   对于连续的两点P1,P2，计算法向量P1t=cross(P1,B), P2t=cross(P2,B)
%%   相位延迟 = arccos(dot(P1t,P2t))
%% ============================================================================
function [Subphr_bg_rm_Ps] = ddgPsPhr(ps_bg_rm,B_bg_rm)
    Avnum = size(ps_bg_rm,1);
    for k = 1:Avnum-1
        P1t_bg_rm = cross(ps_bg_rm(k,:),B_bg_rm);
        P2t_bg_rm = cross(ps_bg_rm(k+1,:),B_bg_rm);
        Subphr_bg_rm_Ps(k) = dot(P1t_bg_rm,P2t_bg_rm)/(sqrt(sum(P1t_bg_rm.^2)*sum(P2t_bg_rm.^2)));
    end
    Subphr_bg_rm_Ps(Subphr_bg_rm_Ps > 1) = 1;
    Subphr_bg_rm_Ps(Subphr_bg_rm_Ps < -1) = -1;
end

%% ============================================================================
%% 子函数: 双折射轴滤波
%% ============================================================================
%% 函数功能:
%%   对3D双折射轴数据进行空间滤波
%%   分离处理Q、U、V三个分量，然后重新归一化
%%
%% 输入参数:
%%   loaxis2: 输入的双折射轴数据 (深度×宽度×3)
%%   h2: 2D滤波器核
%%
%% 输出参数:
%%   loaxis3: 滤波后的双折射轴数据 (深度×宽度×3)
%%
%% 处理步骤:
%%   1. 分离Q、U、V分量
%%   2. 对每个分量应用imfilter滤波
%%   3. 重新归一化向量长度
%% ============================================================================
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

    Temp_ax1(isnan(Temp_ax1)) = 0;
    Temp_ax2(isnan(Temp_ax2)) = 0;
    Temp_ax3(isnan(Temp_ax3)) = 0;
    loaxis3=cat(3,Temp_ax1,Temp_ax2,Temp_ax3);
end

%% ============================================================================
%% 子函数: 深度分辨率恢复 (Depth Resolved LA)
%% ============================================================================
%% 函数功能:
%%   通过连续旋转恢复每个深度位置的双折射轴
%%   基于累积双折射轴和相位延迟进行深度方向的旋转累积
%%
%% 输入参数:
%%   loaxis22: 累积双折射轴 (深度×宽度×3)
%%   phR_gF: 相位延迟 (深度×宽度)
%%   rotationVector2: 旋转向量 (深度×宽度×3)
%%
%% 输出参数:
%%   axiss2: 深度分辨的双折射轴 (深度×宽度×3)
%%
%% 计算原理:
%%   从表面开始，逐深度层累积旋转：
%%   每个深度i的光轴 = 旋转矩阵(i-1) * 光轴(i-1)
%%   旋转角度 = -phR_gF(i)/2
%% ============================================================================
function [axiss2] = drLA(loaxis22,phR_gF,rotationVector2)
    for j=1:size(loaxis22,2)
        for i=1:size(loaxis22,1)
            vector1=squeeze(rotationVector2(1,j,:));
            a=[loaxis22(i,j,1),loaxis22(i,j,2),loaxis22(i,j,3)];
            if i==1
                rotationMatrix1=rotationVectorToMatrix(vector1);
                axis2(i,j,1:3)=a;
            else
                axis2(i,j,1:3)=rotationMatrix1*a';
                Qa=axis2(i,j,1); Ua=axis2(i,j,2); Va=axis2(i,j,3);
                d1= -phR_gF(i,j)/2.0*[Qa,Ua,Va];
                rotationMatrix1 = rotationVectorToMatrix(d1)*rotationMatrix1;
                Q2=Qa/sqrt(Qa^2+Ua^2+Va^2);
                U2=Ua/sqrt(Qa^2+Ua^2+Va^2);
                V2=Va/sqrt(Qa^2+Ua^2+Va^2);
                V2(isnan(V2)) = 0;
                axiss2(i,j,:)=cat(3,Q2,U2,V2);
            end
        end
    end
end
