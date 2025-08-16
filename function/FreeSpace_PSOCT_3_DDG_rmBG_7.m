% 高斯滤波代替中值和平滑。
% QUV小尺度滤波（h1)，光轴大尺度滤波(h2)。
% 对每个A-scan去除DC分量，并且执行标准化
% 用滤波后的光轴计算旋转矩阵。
% output cumulative LA as well ---20240825
% use DDG to cal La and phr
% use T3 as a marker to diff the BG and polar
% correct the BG effect use norm(P3-P1) as a 
% output the Phr with BG removal
% cal dr LA need the Phr with BG induced phr---20240831
% test for retina imaging --------------------20240902
% use more than 3 data point to cal LA and then avg phr
% for volume data

function [LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = FreeSpace_PSOCT_3_DDG_rmBG_7(Qm,Um,Vm,Stru,test_seg_top,h1,h2,Avnum)
showimg = 0;


%% flatting
StruF=Stru*0;QmF=Qm*0;UmF=Um*0;VmF=Vm*0;
for j=1:size(Stru,2)
   StruF(1:size(Stru,1)-test_seg_top(j)+1,j)=Stru(test_seg_top(j):size(Stru,1),j);
   QmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Qm(test_seg_top(j):size(Stru,1),j);
   UmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Um(test_seg_top(j):size(Stru,1),j);
   VmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Vm(test_seg_top(j):size(Stru,1),j);
end


%%  QUV filtering

% QmF=imfilter(QmF,h1,'replicate'); % 个人认为不应该在这里滤波
% UmF=imfilter(UmF,h1,'replicate'); 
% VmF=imfilter(VmF,h1,'replicate');


% QF=QmF-mean(QmF); % zero-meaning
% UF=UmF-mean(UmF);
% VF=VmF-mean(VmF);
QF=QmF;UF=UmF;VF=VmF;

QF1=QF./sqrt(QF.^2+UF.^2+VF.^2); % normallization
UF1=UF./sqrt(QF.^2+UF.^2+VF.^2);
VF1=VF./sqrt(QF.^2+UF.^2+VF.^2);

QF1(isnan(QF1)) = 0;    UF1(isnan(UF1)) = 0;    VF1(isnan(VF1)) = 0;

SmapF=cat(3,QF1,UF1,VF1);
% figure(11),imagesc(SmapF);
%% rot 90 degree
% vetorq=pi/2*[0;1;0];
% rotq=rotationVectorToMatrix(vetorq);
% for i=1:size(SmapF,1)
%     for j=1:size(SmapF,2)
%         aa=squeeze(SmapF(i,j,:));
%         SmapF1(i,j,1:3)=rotq*aa;
%     end
% end

SmapF1=smoothdata(SmapF, 1, 'movmean', 15);
%% svd to fit LA 
stLayer =1;
DTheta = 1.0675*2*pi; %[rad] estimated maximum phr around Q from BG(50 cm more)
dTheta = DTheta/320; %[rad/px]
LABG = [1 0 0];
% Avnum = 3;
phd=zeros(size(SmapF1,1),size(SmapF1,2));
for j=1:size(SmapF1,2)
    dirCheck = 1;
    for i=1:size(SmapF1,1)-Avnum
        planeData=squeeze(SmapF1(i:i+Avnum,j,:)); % QUV points from Avnum layers
%         planeData=squeeze(SmapF1([i,i+Astep,i+Astep*2],j,:)); % QUV points from Avnum layers
        x1=squeeze(planeData(:,1)); y1=squeeze(planeData(:,2)); z1=squeeze(planeData(:,3));
        nx=sum(sum(x1==0));
        if x1==0|nx>0 % skip if Qs has 0
            cumLA_bg(i,j,1:3)=0;
            cumLA_raw(i,j,1:3)=0;
            phR(i,j)=0;
            phR_rmBG(i,j) = 0;
            phR_raw(i,j) = 0;
            phd(i,j)=0;
            rotationVector(i,j,1:3)=0;
        else
            % rotate the det p to remove DC part of BG
            ps_bg_rdu = [];ps_bg_rm = [];
            for jL =1:Avnum
                inpp = planeData(jL,:);
                outp = rotationVectorToMatrix(-1*(i+test_seg_top(j)+stLayer-0)*dTheta*LABG)*inpp';
                outp2 = rotationVectorToMatrix(-1*(i+test_seg_top(j)+stLayer+jL-2)*dTheta*LABG)*inpp';
%                 outp = rotationVectorToMatrix((i+test_seg_top(j)+stLayer-0)*dTheta*LABG)*inpp';
%                 outp2 = rotationVectorToMatrix((i+test_seg_top(j)+stLayer+jL-2)*dTheta*LABG)*inpp';
                ps_bg_rdu = [ps_bg_rdu;outp'];
                ps_bg_rm = [ps_bg_rm;outp2'];
            end
            
            % cal the LA while correct the BG: prepare points and vectors
            % [DDG: pick points at start, middle and end to cal LA]
            P1  = ps_bg_rm(1,:); P2 = ps_bg_rm(round((1+Avnum)/2),:); P3 = ps_bg_rm(Avnum,:);
            P1p = ps_bg_rdu(1,:);P2p = ps_bg_rdu(2,:);P3p = ps_bg_rdu(3,:);
            T1p = P2p - P1p; T2p = P3p-P2p; T1 = P2-P1;T2 = P3-P2;
            Tb1 = P1 -P1p; Tb2 = P3-P3p;
            
            % direct cal LA and phr from bg rm points [bgrm DDG]
            B_bg_rm = cross(T1,T2); B_bg_rm = B_bg_rm/sqrt(sum(B_bg_rm.^2));
            % ddg use normal vector of T to cal phr between adjacent points
            [Subphr_bg_rm_Ts] = ddgTsPhr(ps_bg_rm,B_bg_rm);
% %             phr_bg_rm_Ts(i,j) = median(acos(Subphr_bg_rm_Ts));% acos first , then avg
%             phr_bg_rm_Ts(i,j) = acos(median(Subphr_bg_rm_Ts));% avg acos first, then acos
            
            % ddg use tangent vector to cal phr between adjacent points
%             [Subphr_bg_rm_Ps] = ddgPsPhr(ps_bg_rm,B_bg_rm);
%             phR_bgrm_Ps2(i,j) = median(acos(Subphr_bg_rm_Ps));% acos first , then avg
            phR_raw(i,j) = acos(median(Subphr_bg_rm_Ts));% avg acos first, then acos
            cumLA_raw(i,j,1:3) = -B_bg_rm;
            
            % calc Bp with bg and then correct Bp with deltaB [BG corr DDG]
            Bp = cross(T1p,T2p); 
            T3_nor = norm(P3-P1);
            % condition to check the BG and cal deltaB
            if T3_nor < 0.065%0.7/0.12
                deltaB = [0 0 0];
                T1 = T1p; T2 = T2p;
                if dirCheck && (Bp(1) > 0) %% check LA direction
                    Bp = -Bp;
                end
%                 if i == 10, Bp = -Bp;end %% for 2024.07.03_15.54.54_2.zhaoyu_repB only
            else
                dirCheck = 0;
                deltaB = cross(T1p,Tb2) + cross(T2p,Tb1) - cross(Tb1,Tb2);% BG correction
            end
            B = Bp + deltaB; % corrected LA
            B = B/sqrt(sum(B.^2)); % normalize
            T3_nors(i,j) = T3_nor;
            if 1
                if T3_nor < 0.065
                    % cal phr use Ps ()
%                     phR(i,j) = acos(median(ddgPsPhr(ps_bg_rdu,B)));
                    phR(i,j) = acos(median(ddgTsPhr(ps_bg_rdu,B)));
                    phR_rmBG(i,j) = max([phR(i,j) - 1*dTheta, 0]);
                    phR_rmBG(i,j) = min([phR(i,j), 0.2]);
                else
%                     phR(i,j) = acos(median(ddgPsPhr(ps_bg_rm,B)));
                    phR(i,j) = acos(median(ddgTsPhr(ps_bg_rm,B)));
                    phR_rmBG(i,j) = phR(i,j);
                end
                cumLA_bg(i,j,1:3)=-B;
            else
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
% post process the LA 
%% 光轴滤波
[cumLA_bg_gF] = LAgFfilt(cumLA_bg,h2);
[cumLA_raw_gF] = LAgFfilt(cumLA_raw,h2);
if showimg
figure;imshow(cumLA_bg_gF*0.5+0.5,[])
figure;imshow(cumLA_raw_gF*0.5+0.5,[])
figure;imshow(phR_rmBG,[0 0.25])
figure;imshow(phR_raw,[0 0.25])
end
%% 相位滤波
    phR_gF=imfilter(phR,h2,'replicate');%% with BG
    phR_rmBG_gF=imfilter(phR_rmBG,h2,'replicate');%% with BG
    phR_raw_gF = imfilter(phR_raw,h2,'replicate');%% raw w/o BG
    phrRg = [0 0.3];
    if showimg
    figure;imshowpair(mat2gray(phR_rmBG,phrRg),mat2gray(phR_rmBG_gF,phrRg),'montage')
     figure;imshowpair(mat2gray(phR_raw,phrRg),mat2gray(phR_raw_gF,phrRg),'montage')
    end
%% 计算旋转矩阵
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


%% rotate the LA to achieve the depth resolved L
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
% restore curved results
LA=bfloaxis3D*0;PhR=bfphrr*0;cumLA = bfloaxis3D*0;
LA_raw=bfloaxis3D*0;PhR_raw=bfphrr*0;cumLA_raw = bfloaxis3D*0;
for j=1:size(Stru,2)
    toInd = test_seg_top(j):size(LA,1); fromInd = 1:size(LA,1)-test_seg_top(j)+1;
    LA(toInd,j,:)=bfloaxis3D(fromInd,j,:);
    PhR(toInd,j)=bfphrr(fromInd,j);
    cumLA(toInd,j,:)=cumLA_bg_gF(fromInd,j,:);
    
    LA_raw(toInd,j,:)=bfloaxis3D_raw(fromInd,j,:);
    PhR_raw(toInd,j)=bfphrr_raw(fromInd,j);
    cumLA_raw(toInd,j,:)=cumLA_raw_gF(fromInd,j,:); 
end
end

%% cal phr used Ts: P2P1 and P3P2
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
%% cal phr used Ps: P1 and P2
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

%% filt cumLA
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
%% rot to achieve drLA
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
