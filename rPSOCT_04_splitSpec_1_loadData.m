%% ini params: streamline forLoop from rPSOCT_03_splitSpec_2_repB
clear all;close all;
clc;
if ~exist('.\output\','dir'),mkdir('.\output\'),end
outputdir=[cd,'\output\'];
disp_coef=-20.1;
do_PhComp=1; do_medianshift=1;
do_reg=0; % toggle Inter-Frame registration (for motion removal)
useref=1;
saveDicom=0;
show_img=1;


distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing
showimg=1;
tmp=dir('.\*.oct');
for iFile=1:size(tmp,1), fprintf(['[%d]\t',tmp(iFile).name(1:end-4),'\n'],iFile);end
%% choose file and select Z range
tic
iZShift = zeros(1,numel(tmp));
iy=1;
do_ed = 0;

%% 参数设置
Avnum = 7;   %number of layers made use in first calculation to compute the normal vector in Poincare sphere
h1 = fspecial('gaussian',[7 7],5);%预处理
h2 = fspecial('gaussian',[21 21],9);%预处理

for iFile=1%
    Thr=170; % Threshold to remove low OCT signal
    % This reads the parameter used for the data acquisition from *.oct* file
    filename=tmp(iFile).name;
    fid=fopen(filename);
    bob=fread(fid,1,'uint32');
    SPL=fread(fid,1,'double');
    nX=fread(fid,1,'uint32'); %% number of Alines
    nY=fread(fid,1,'uint32'); %% number of B-scans
    Boffset=fread(fid,1,'uint32');
    Blength=fread(fid,1,'uint32')+1;
    Xcenter=fread(fid,1,'double');
    Xspan=fread(fid,1,'double'); 
    Ycenter=fread(fid,1,'double');
    Yspan=fread(fid,1,'double');
    frame_per_pos=fread(fid,1,'uint32'); %% repetition of B-scans
    n_dataset=fread(fid,1,'uint32'); %% repetion of volume scan
    ProtMode=fread(fid,1,'uint32');
    fseek(fid,4,'cof');%v10
    sizeBck=fread(fid,1,'uint32');
    Bck1=fread(fid,sizeBck,'int16');
    sizeKES=fread(fid,2,'uint32');
    KES1=(fread(fid,sizeKES(2),'double'))'*sizeKES(2);
    sizeBck=fread(fid,1,'uint32');
    Bck2=fread(fid,sizeBck,'int16');
    sizeKES=fread(fid,2,'uint32'); %new 4096 cal
    KES2=(fread(fid,sizeKES(2),'double'))'*sizeKES(2);
    disp_coef=fread(fid,1,'double'); %%dispersion coefficient
    nR=frame_per_pos;
    IMGheight=floor(Blength/2);
    kmat=linspace(-0.5,0.5,Blength)'.^2;
    phV=exp(1i.*(kmat.*disp_coef));
    
    nY=floor(nY/nR);
        % looks likes these param unsed
    K=1:Blength;
    use_autoRg=1;
    RgFlow=[60 110];
    jsurf=zeros(4);
    IMcropRg=1:(SPL/2);
    nZcrop=numel(IMcropRg);
    imshowrgZ=1:nZcrop;
    jusamp=zeros(nZcrop,nX,4);
    winG =   tukeywin(Blength,0.25);
    % subWins 鍙傛暟
    nWin = 9; 
    winL = 2*Blength/(nWin+1);
    winG=tukeywin(winL,0.25);
    winG_whole = tukeywin(Blength,0.25); % window for whole spectrum
    windex=1:winL/2:Blength;
    
    if useref==1 %use the first 50k a-lines to calc ref
        fseek(fid,bob,'bof');
        n_ref=floor(min(50000,floor(nY*nX/2)*2)/nX);
        Ref_ch1=zeros(SPL,nX,n_ref);
        Ref_ch2=zeros(SPL,nX,n_ref);
        for i_ref=1:n_ref
            fseek(fid,4 ,'cof');
            BB=reshape(fread(fid,SPL*nX*2,'int16'),[SPL,nX*2]);
            Ref_ch1(:,:,i_ref)=BB(:,1:nX);
            Ref_ch2(:,:,i_ref)=BB(:,nX+1:end);
        end
        Ref_ch1=repmat(mean(mean(Ref_ch1,2),3),[1,nX]);
        Ref_ch2=repmat(mean(mean(Ref_ch2,2),3),[1,nX]);
    elseif useref==-1
        Ref_ch1=0;
        Ref_ch2=0;
    else
        Ref_ch1=repmat(Bck1,[1,nX]);
        Ref_ch2=repmat(Bck2,[1,nX]);
    end    

   
%% Calculates the Cumulative Stokes parameters I,Q,U,V 
for nr = [nR]
foutputdir = fullfile(outputdir,[filename,'.rep',num2str(nr),'\']);
mkdir(foutputdir); 
czrg = [51:500]+iZShift(iFile);
czrg = [1:320];
Qms = zeros(numel(czrg),nX,nY);
Ums = Qms;Vms = Qms;Strus = Qms;
SSmap = zeros(numel(czrg),nX,3,nY);
pb = ProgressBar(nY);
dopu_splitSpec = zeros(nZcrop,nX,nY,nr);
% parfor iY=1:iy:nY % Bscan number
%     for iY = 1:iy:nY
for iY = 260 % ===============================>>>>>>>>>>> 鎸囧畾 绗? 100 甯?
        
        Bs1=zeros(Blength,nX,nr);
        Bs2=zeros(Blength,nX,nr);
        fid=fopen(filename);
        fseek(fid,bob+(SPL*nX*2+2)*2*(nR*(iY-1)),'bof');
        for Ic=1:nr
            fseek(fid,4,'cof');  % 4byte =2 int16
            B=double(reshape(fread(fid,SPL*nX*2,'int16'),[SPL,nX*2]));
            Bs1(:,:,Ic) = (B(:,1:nX) - Ref_ch1);
            Bs2(:,:,Ic) = (B(:,nX+1:end) - Ref_ch2);
        end
        if do_PhComp==1
%             Bd1 =hilbert(Bs1).*phV;
%             Bd2 =hilbert(Bs2).*phV;
            Bd1 =real(hilbert(Bs1).*phV);  % dispersion correction
            Bd2 =real(hilbert(Bs2).*phV);
        else
            Bd1=Bs1;
            Bd2=Bs2;
        end
        %% split spectrum
        Bimg1 = zeros(SPL,nX,nr,nWin);Bimg2 = Bimg1;
        S0 = zeros(nZcrop,nX,nr,nWin); S1 = S0;S2=S0;S3=S0;
        for iR = 1:nr
            for iL=1:nWin
                iBd1=Bd1(windex(iL):windex(iL)+winL-1,:,iR).*winG;
                iBd2=Bd2(windex(iL):windex(iL)+winL-1,:,iR).*winG;
                Bimg1(:,:,iR,iL)=fft(iBd1,SPL,1);
                Bimg2(:,:,iR,iL)=fft(iBd2,SPL,1);
            end
        end
        IMGs_ch1=Bimg1(IMcropRg,:,:,:);IMGs_ch2=Bimg2(IMcropRg,:,:,:);
        for iR = 1:nr
            for iL = 1: nWin
                [S0(:,:,iR,iL),S1(:,:,iR,iL),S2(:,:,iR,iL),S3(:,:,iR,iL)] = ...
                    cumulativeQUV(IMGs_ch1(:,:,iR,iL),IMGs_ch2(:,:,iR,iL));
            end
        end
        dopu_splitSpec(:,:,iY,:) = sqrt(mean(S1./S0,4).^2+mean(S2./S0,4).^2+mean(S3./S0,4).^2);%% ===>
        dopu_splitSpec_M(:,:)=mean(dopu_splitSpec(:,:,iY,:),4);
        %% whole spectrum nZ*nX*nR
        Bimg1_wholeStr=fft(Bd1.*winG_whole,SPL,1);Bimg2_wholeStr=fft(Bd2.*winG_whole,SPL,1);
        IMG1_wholeStr = Bimg1_wholeStr(IMcropRg,:,:);IMG2_wholeStr = Bimg2_wholeStr(IMcropRg,:,:);
        whlCh1(:,:,iY,:) = IMG1_wholeStr; whlCh2(:,:,iY,:) = IMG2_wholeStr;%% ===>
        %%
%         figure;imshow(20*log10(abs(IMG2_wholeStr(:,:,1))),[80 120])
        fclose(fid);
        pb.progress;
end
pb.stop;

fclose all;
end
end
%% ==>>>output: whlCh1(:,:,iY,:),whlCh2,dopu_splitSpec(:,:,iY,iR,iL)<<<<==============================

%% output variable: 
% whlCh1: Fourier transform of Channel 1
% whlCh2: Fourier transform of Channel 2

%% take 1 frame from whlCh1 and whlCh2
% iY = 100; 
ch1 = squeeze(whlCh1(:,:,iY,:));
ch2 = squeeze(whlCh2(:,:,iY,:));
Stru=20*log10(abs(ch1(:,:,1)).^2+abs(ch2(:,:,1)).^2);
figure(1);imshow(Stru,[170 250])

Stru_orig=abs(ch1(:,:,1)).^2+abs(ch2(:,:,1)).^2;
Stru_M=Stru_orig;   % 变核加权滤波

%% OAC计算
Nois_M=mean(mean(Stru_M(size(Stru_M,1)-2:size(Stru_M,1),:)));
Nois_D=std(std(Stru_M(size(Stru_M,1)-2:size(Stru_M,1),:)));
Stru_M=Stru_M-Nois_M+Nois_D; %去除噪声地板
down=Stru_M<0;Stru_M(down)=1; %去除负数部分
tail_signal=mean(Stru_M(size(Stru_M,1)-5:size(Stru_M,1),:)); %估计剩余光强
tail_signal = medfilt1(tail_signal,15);
tail_signal = smooth(tail_signal, 15)';

OAC=Stru_M*0;
for z=1:size(Stru_M,1)
    OAC(z,:)=Stru_M(z,:)./(2*0.005*sum(Stru_M(z+1:size(Stru_M,1),:))+tail_signal);
end
% figure(2),imshow(OAC,[0 3]);

%% 表面检测(如果有label可以直接载入)

test_seg_top=surf_seg(OAC,0.1);
figure(3);imshow(Stru,[170 250]);hold on 
plot(test_seg_top,'-r', 'LineWidth', 2);


%% eigen to extract flow
czrg = 1:size(ch1,1); % depth range
[p_tis,p_bld] = OCTA_F_ED_Clutter_EigFeed(cat(1,ch1,ch2), 2);%提取动态成分
p_tis1 = 20*log10(p_tis(1:numel(czrg),:));
p_tis2 = 20*log10(p_tis(numel(czrg)+1:end,:));
p_bld1 = 20*log10(p_bld(1:numel(czrg),:));
p_bld2 = 20*log10(p_bld(numel(czrg)+1:end,:));    
BL1 = p_bld1 - 0.5*imgaussfilt(p_tis1,5);
BL2 = p_bld2 - 0.5*imgaussfilt(p_tis2,5);
BL = BL1+BL2;
figure(4);imshow(BL,[70 120])
%% Average 
kRL=2;kRU=7;
A_ch1=mean(ch1,3);
A_ch2=mean(ch2,3);
[AS0,AS1,AS2,AS3] = cumulativeQUV(A_ch1,A_ch2);
AQm=AS1./AS0;
AUm=AS2./AS0;
AVm=AS3./AS0;
QUV_A(:,:,:)=cat(3,AQm,AUm,AVm);
Stru_A=20*log10(S0);

[AQmm] = vWinAvgFiltOpt_2_1(AQm,dopu_splitSpec_M,kRL,kRU);
[AUmm] = vWinAvgFiltOpt_2_1(AUm,dopu_splitSpec_M,kRL,kRU);
[AVmm] = vWinAvgFiltOpt_2_1(AVm,dopu_splitSpec_M,kRL,kRU);

[LA_A,PhR_A] = FreeSpace_PSOCT(AQmm,AUmm,AVmm,Stru_A,test_seg_top,h1,h2,Avnum);
figure(5),imagesc((QUV_A+1)/2);
figure(6),imagesc((LA_A+1)/2);
figure(7),imagesc(PhR_A,[0 0.5]);

%% eigen to extract static components
[IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,ch1,ch2), 1); % 提取静态分量?
IMG_ch1 = IMG_ch(1:size(ch1,1),:);IMG_ch2 = IMG_ch(size(ch1,1)+1:end,:);
%figure(2);imshow(20*log10(abs(IMG_ch1).^2+abs(IMG_ch2).^2),[170 250])

[ES0,ES1,ES2,ES3] = cumulativeQUV(IMG_ch1,IMG_ch2);
EQm=ES1./ES0;
EUm=ES2./ES0;
EVm=ES3./ES0;
QUV_E(:,:,:)=cat(3,EQm,EUm,EVm);
Stru_E=20*log10(S0);

[EQmm] = vWinAvgFiltOpt_2_1(EQm,dopu_splitSpec_M,kRL,kRU);
[EUmm] = vWinAvgFiltOpt_2_1(EUm,dopu_splitSpec_M,kRL,kRU);
[EVmm] = vWinAvgFiltOpt_2_1(EVm,dopu_splitSpec_M,kRL,kRU);

[LA_E,PhR_E] = FreeSpace_PSOCT(EQmm,EUmm,EVmm,Stru_E,test_seg_top,h1,h2,Avnum);
figure(8),imagesc((QUV_E+1)/2);
figure(9),imagesc((LA_E+1)/2);
figure(10),imagesc(PhR_E,[0 0.5]);

% figure(11),imshow(AQm./EQm,[0.7 1.3]);
% figure(12),imshow(AUm./EUm,[0.7 1.3]);
% figure(13),imshow(AVm./EVm,[0.7 1.3]);
figure(14),imshow(dopu_splitSpec_M,[]);
%%
% cal QUV
function [S0,S1,S2,S3] = cumulativeQUV(IMG_ch1,IMG_ch2)
axis=angle(IMG_ch2.*conj(IMG_ch1));
S0=abs(IMG_ch1).^2+abs(IMG_ch2).^2;
S1=abs(IMG_ch1).^2-abs(IMG_ch2).^2;
S2=2.*abs(IMG_ch1).*abs(IMG_ch2).*cos(axis);
S3=2.*abs(IMG_ch1).*abs(IMG_ch2).*sin(-axis);
end


