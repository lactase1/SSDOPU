%% plotFigures: figure5-2
%% ini params
close all;clear;clc;
DTheta = 1.0675*2*pi; %[rad] estimated maximum phr around Q from BG(50 cm more)
dTheta = DTheta/320; %[rad/px] 0.23deg/um
LA1 = [-1, 1, -1];
LA1 = LA1./sqrt(sum(LA1.^2))
figure;
[X Y Z] = sphere(20);
SS=surf(X,Y,Z);
SS.FaceAlpha = 0.05;
SS.EdgeAlpha = 0.05;
xlim([-1.1 1.1])
axis equal
xlabel('Q')
ylabel('U')
zlabel('V')
grid on;hold on
plot3([-LA1(1),LA1(1)],[-LA1(2),LA1(2)],[-LA1(3),LA1(3)],'or-.','linewidth',3)
%% rot around LA1 [polarized material]
Np = 30;
Np_curve = 100;
phr1 = 0.07;%[rad/px] 0.8 deg/um *5*3.14/180
p1 = [-0.8,-0.4,0.24];
p1 = p1/sqrt(sum(p1.^2));
p_slab1 = p1;
for i = 1:Np
    inpp = p_slab1(i,:);
    outp = rotationVectorToMatrix(phr1*LA1)*inpp';
    p_slab1 = [p_slab1;outp'];
end

p_slab1_nonP = repmat(p1,Np+1,1);%% nonPolarized material [w/o noise]
plot3(p_slab1(:,1),p_slab1(:,2),p_slab1(:,3),'ok','linewidth',1,'MarkerFaceColor','r','MarkerSize',4)
[p_slab1_curve] = gencirCurve(p1,LA1,Np_curve);
plot3(p_slab1_curve(:,1),p_slab1_curve(:,2),p_slab1_curve(:,3),'r-','linewidth',1)

% return
%% rot around BG
LABG = [1 0 0];
LABG = LABG./sqrt(sum(LABG.^2));
p_slab1_BG = p1;
for i = 2:Np
    inpp = p_slab1(i,:);
    outp = rotationVectorToMatrix((i-1)*dTheta*LABG)*inpp';
    [rot_curve] = gencirCurve(inpp,LABG,Np_curve,(i-1)*dTheta);
    p_slab1_BG = [p_slab1_BG;outp'];
end
[p_slab1_BG] = stokesMixBG(p_slab1,dTheta,LABG);
[p_slab1_nonP_BG] = stokesMixBG(p_slab1_nonP,dTheta,LABG);
plot3([-LABG(1),LABG(1)],[-LABG(2),LABG(2)],[-LABG(3),LABG(3)],'om-.','linewidth',3)
plot3(p_slab1_BG(:,1),p_slab1_BG(:,2),p_slab1_BG(:,3),'ok-','linewidth',1,'MarkerFaceColor','g','MarkerSize',4)
view(-30,0)
%% DDG fitting
avgNum = 3;
for i = [10]
    % ddg cal cum LA and phr of det signal w/o BG
    [LA(i,:),phrs(i)] = cumDDGLaPhr(p_slab1(i:i+avgNum-1,:));
    %rotate det ps back to close the raw points
    ps_bg = p_slab1_BG(i:i+avgNum-1,:);
    ps_bg_rdu = [];ps_bg_rm = [];
    for j =1: avgNum
        inpp = ps_bg(j,:);
        outp = rotationVectorToMatrix(-(i-0)*dTheta*LABG)*inpp';
        outp2 = rotationVectorToMatrix(-(i+j-2)*dTheta*LABG)*inpp';
        [rot_curve] = gencirCurve(inpp,LABG,Np_curve,-(i-0)*dTheta);
        plot3(rot_curve(:,1),rot_curve(:,2),rot_curve(:,3),'m-','linewidth',1)
        ps_bg_rdu = [ps_bg_rdu;outp'];
        ps_bg_rm = [ps_bg_rm;outp2'];
    end
    [ps_bg_rdu,ps_bg_rm] = stokesBGrdm(p_slab1_BG(i:i+avgNum-1,:),i,dTheta,LABG,avgNum);
    plot3(ps_bg_rdu(:,1),ps_bg_rdu(:,2),ps_bg_rdu(:,3),...
        'ok-','linewidth',1,'MarkerFaceColor','b','MarkerSize',4);% reduced BG singal
    %
    [LA_wbg(i,:),~] = cumDDGLaPhr(ps_bg_rdu);
    % ddg cal correctted cum LA and phr from BG
    [B(i,:), phrs_p(i)] = cumDDGLaPhrBG(ps_bg_rm,ps_bg_rdu);
    % nonP w/BG, wo noise
    [ps_bg_rdu,ps_bg_rm] = stokesBGrdm(p_slab1_nonP_BG(i:i+avgNum-1,:),i,dTheta,LABG,avgNum);
    [B_nonp(i,:), phrs_nonp(i)] = cumDDGLaPhrBG(ps_bg_rm,ps_bg_rdu);
    % nonP w/BG, w/ noise
    [B_nonpnoise(i,:), phrs_nonpnoise(i)] = cumDDGLaPhrBG(ps_bg_rm,ps_bg_rdu);
    plot3([-LA_wbg(i,1),LA_wbg(i,1)],[-LA_wbg(i,2),LA_wbg(i,2)],[-LA_wbg(i,3),LA_wbg(i,3)],'b-','linewidth',1.5)
    plot3([-B(i,1),B(i,1)],[-B(i,2),B(i,2)],[-B(i,3),B(i,3)],'k-','linewidth',1.5)
end
%%

% generate a circle of a point around a LA
function [p_slab1_curve] = gencirCurve(p1,LA1,Np_curve,theta)
THrg = 2*pi;
if nargin > 3, THrg = theta;end
p_slab1_curve = p1;
for i = 1:Np_curve
    inpp = p_slab1_curve(i,:);
    outp = rotationVectorToMatrix(THrg/Np_curve*LA1)*inpp';
    p_slab1_curve = [p_slab1_curve;outp'];
end
end

% DDG to cal LA and phr
function [LA,phr] = cumDDGLaPhr(ps3)
P1 = ps3(1,:);P2 = ps3(2,:);P3 = ps3(3,:);
T1 = (P2-P1); T1 = T1/sqrt(sum(T1.^2));
T2 = (P3-P2); T2 = T2/sqrt(sum(T2.^2));
B1 = cross(T1,T2); B1 = B1/sqrt(sum(B1.^2));
N1 = cross(B1,T1);N2 = cross(B1,T2);
delta = 2*0.5*acos(dot(N1,N2)/(sqrt(sum(N1.^2)*sum(N2.^2))));
LA = B1; phr = delta;
end
% DDG to cal LA and phr with BG correction
function [LA, phr] = cumDDGLaPhrBG(ps_bg_rm,ps_bg_rdu)
P1  = ps_bg_rm(1,:); P2 = ps_bg_rm(2,:); P3 = ps_bg_rm(3,:);
P1p = ps_bg_rdu(1,:);P2p = ps_bg_rdu(2,:);P3p = ps_bg_rdu(3,:);
T1p = P2p - P1p; T2p = P3p-P2p; T1 = P2-P1;T2 = P3-P2;
Tb1 = P1 -P1p; Tb2 = P3-P3p;
Bp = cross(T1p,T2p);
T3_nor = norm(P3-P1);
if T3_nor < 0.12
    deltaB = [0 0 0];
    T1 = T1p; T2 = T2p;
else
    deltaB = cross(T1p,Tb2) + cross(T2p,Tb1) - cross(Tb1,Tb2);% BG correction
end
B = Bp + deltaB; % corrected LA
B = B/sqrt(sum(B.^2)); % normalize
N1 = cross(B,T1); N2 = cross(B,T2);% normal vector
phr = acos(dot(N1,N2)/sqrt(sum(N1.^2)*sum(N2.^2)));% phase retardation
LA = -B;
end
% recover the BG effect on Stokes
function [ps_bg_rdu,ps_bg_rm] = stokesBGrdm(ps_bg,iLayer,dTheta,LABG,avgNum)
ps_bg_rdu = [];ps_bg_rm = [];
for j =1: avgNum
    inpp = ps_bg(j,:);
    outp = rotationVectorToMatrix(-(iLayer-0)*dTheta*LABG)*inpp';
    outp2 = rotationVectorToMatrix(-(iLayer+j-2)*dTheta*LABG)*inpp';
    ps_bg_rdu = [ps_bg_rdu;outp'];
    ps_bg_rm = [ps_bg_rm;outp2'];
end
end
% stokes mixed BG
function [p_slab1_BG] = stokesMixBG(p_slab1,dTheta,LABG)
p_slab1_BG = p_slab1(1,:);
for i = 2:size(p_slab1,1)
    inpp = p_slab1(i,:);
    outp = rotationVectorToMatrix((i-1)*dTheta*LABG)*inpp';
    p_slab1_BG = [p_slab1_BG;outp'];
end
end

% 增加输出路径控制
outputPath = fullfile(pwd, 'output'); % 默认输出路径为当前目录下的 output 文件夹
if ~exist(outputPath, 'dir')
    mkdir(outputPath); % 如果路径不存在，则创建
end

% 保存图形
saveas(gcf, fullfile(outputPath, 'figure5_2.png')); % 保存为 PNG 格式
savefig(gcf, fullfile(outputPath, 'figure5_2.fig')); % 保存为 MATLAB FIG 格式
