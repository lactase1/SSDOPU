function rPSOCT_process_single_file_gpu(varargin)
% This function processes a single rPSOCT data file using GPU acceleration if available.
% 
% Usage:
%   rPSOCT_process_single_file_gpu('ParameterName', ParameterValue, ...)
%
% GPU-accelerated version based on rPSOCT_06_splitSpec_3_3D_without_mat_wyx.m
%
% Author: [Yongxin Wang]
% Date: [2023]

%% 检测是否可以使用GPU
try
    gpuDevice();
    fprintf('GPU可用，将使用GPU加速计算\n');
    gpu_available = true;
catch
    fprintf('GPU不可用或MATLAB未安装Parallel Computing Toolbox，将使用CPU计算\n');
    gpu_available = false;
end

%% 解析输入参数
p = inputParser;
% 添加参数
addParameter(p,'file','',@ischar);
addParameter(p,'outDir','',@ischar);
addParameter(p,'display_name','',@ischar);
addParameter(p,'do_avg',true,@islogical);
addParameter(p,'do_eig',false,@islogical);
addParameter(p,'do_cfg1',true,@islogical);
addParameter(p,'do_cfg2',false,@islogical);
addParameter(p,'do_PhComp',true,@islogical);
addParameter(p,'do_ssdopu',true,@islogical);
addParameter(p,'show_img',false,@islogical);
addParameter(p,'saveDicom',true,@islogical);
% 解析输入
parse(p,varargin{:});
% 提取参数
f = p.Results.file;
outDir = p.Results.outDir;
display_name = p.Results.display_name;
do_avg = p.Results.do_avg;
do_eig = p.Results.do_eig;
do_cfg1 = p.Results.do_cfg1;
do_cfg2 = p.Results.do_cfg2;
do_PhComp = p.Results.do_PhComp;
do_ssdopu = p.Results.do_ssdopu;
show_img = p.Results.show_img;
saveDicom = p.Results.saveDicom;

%% 检查文件有效性
if ~exist(f,'file')
    error('输入文件 %s 不存在', f);
end

% 创建输出目录
if ~exist(outDir,'dir')
    mkdir(outDir);
end
foutputdir = [outDir,filesep];

% 获取文件名
[~, name, ~] = fileparts(f);
if isempty(display_name)
    display_name = name;
end
fprintf('处理文件: %s\n', display_name);

%% 加载元数据
try
    filename = f;
    fid=fopen(filename);
    seg=fread(fid,4,'int16');
    bob=ftell(fid);
    Blength = double(seg(1)); nR = double(seg(2)); 
    nX = double(seg(3)); nY = double(seg(4));
    fclose(fid);
catch ME
    fprintf('文件元数据读取错误: %s\n', ME.message);
    return;
end

%% 设置各种处理参数
hasSeg = false;
load('DispComp.mat'); % 加载相位校正
winLen = 300; % 窗口长度
AveragingNum = 5; % averaging numer along x

SPL = Blength; 
nr = nR;   % 所需B扫描的数量
Avnum = AveragingNum; % averaging number along z
avgSmethod = 1; % avearge the spectra instead of the processed OA,Ret,biref  1:yes 0:no

wovWinF = true; % windowless variable-window filtering is Ture;otherwise is FALSE
kRL_cfg1 = 7; kRU_cfg1 = 70; % config1
kRL_cfg2 = 7; kRU_cfg2 = 28; % config2
h1 = 3; h2 = 8; % config1/2

%% 加载参考光
Ref_ch1 = double(mean(importdata('Ref_ch1.mat')));
Ref_ch2 = double(mean(importdata('Ref_ch2.mat')));

%% 设置频谱分割参数
nWin = 5;
winL = 700;
windex =[ 1 301 601 901 1201];
g = hann(winL);
winG=repmat(g,[1 nX]);
g_whole = hann(SPL);
winG_whole=repmat(g_whole,[1 nX]);

%% 数据处理
czrg = 2:401; % 截取的深度范围

for nr = [nR]
    % 设置一些相关参数
    nZcrop = numel(czrg);nZ = nZcrop;
    
    %% 为了简化代码，改用直接循环处理
    % 以下代码不再使用parfor或cell数组方式，直接使用全局数组处理

    % 为最终结果预分配数组
    Strus = zeros(numel(czrg),nX,nY);
    Stru_OAC = zeros(numel(czrg),nX,nY);
    dopu_splitSpectrum = zeros(numel(czrg),nX,nY);
    Smap_avg = zeros(numel(czrg),nX,3,nY);
    Smap_rep1 = zeros(numel(czrg),nX,3,nY);

    LA_c_cfg1_avg = zeros(nZcrop-Avnum,nX,3,nY);
    PhR_c_cfg1_avg = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg1_avg = zeros(nZcrop-Avnum,nX,3,nY);
    LA_c_cfg1_eig = zeros(nZcrop-Avnum,nX,3,nY);
    PhR_c_cfg1_eig = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg1_eig = zeros(nZcrop-Avnum,nX,3,nY);

    LA_Ms_cfg1_rmBG = zeros(nZcrop-Avnum,nX,3,nY);
    PhR_Ms_cfg1_rmBG = zeros(nZcrop-Avnum,nX,nY);
    cumLA_Ms_cfg1_rmBG = zeros(nZcrop-Avnum,nX,3,nY);

    LA_c_cfg2_avg = zeros(nZcrop-Avnum,nX,3,nY);
    PhR_c_cfg2_avg = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg2_avg = zeros(nZcrop-Avnum,nX,3,nY);
    LA_c_cfg2_eig = zeros(nZcrop-Avnum,nX,3,nY);
    PhR_c_cfg2_eig = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg2_eig = zeros(nZcrop-Avnum,nX,3,nY);

    pb = ProgressBar(nY);

    % 如果GPU可用，预加载GPU函数
    if gpu_available
        cumulativeQUV_gpu(complex(gpuArray([1 1]), gpuArray([1 1])), complex(gpuArray([1 1]), gpuArray([1 1])));
        vWinAvgFiltOpt_2_1_gpu(gpuArray([1 1]), gpuArray([1 1]), 2, 7);
    end

    % 预处理GPU版本的窗口函数
    gpu_winG = winG;
    gpu_winG_whole = winG_whole;
    if gpu_available
        gpu_winG = gpuArray(winG);
        gpu_winG_whole = gpuArray(winG_whole);
    end

    % 为表面检测预分配内存
    topLines = zeros(nX, nY);

    for iY=1:nY-1 % Bscan number
        % 读取B扫描数据
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
        
        % 如果GPU可用，将数据转移到GPU
        if gpu_available
            Bs1 = gpuArray(Bs1);
            Bs2 = gpuArray(Bs2);
        end
        
        if do_PhComp==1
            Bd1 =real(hilbert(Bs1).*phV);  % dispersion correction
            Bd2 =real(hilbert(Bs2).*phV);
        else
            Bd1=Bs1;
            Bd2=Bs2;
        end
        
        %% split spectrum
        if ~do_ssdopu % if not do adaptive gaussian filter dopu_ss =1 
            dopu_ss = 1;
        else
            % creat array to store split-spectrum complex (FFT): Z*X*nR*nWin 
            Bimg1 = zeros(SPL,nX,nr,nWin);Bimg2 = Bimg1;
            S0 = zeros(nZcrop,nX,nr,nWin); S1 = S0;S2=S0;S3=S0;
            
            % 如果GPU可用，将数组转移到GPU
            if gpu_available
                Bimg1 = gpuArray(Bimg1);
                Bimg2 = gpuArray(Bimg2);
                S0 = gpuArray(S0);
                S1 = gpuArray(S1);
                S2 = gpuArray(S2);
                S3 = gpuArray(S3);
            end
            
            for iR = 1:nr
                for iL=1:nWin
                    % extract data fragments from two different channel and
                    % apply a giassian filter before performing FFT
                    iBd1=Bd1(windex(iL):windex(iL)+winL-1,:,iR).*gpu_winG;
                    iBd2=Bd2(windex(iL):windex(iL)+winL-1,:,iR).*gpu_winG;
                    Bimg1(:,:,iR,iL)=fft(iBd1,SPL,1);
                    Bimg2(:,:,iR,iL)=fft(iBd2,SPL,1);
                end
            end
            IMGs_ch1=Bimg1(czrg,:,:,:);IMGs_ch2=Bimg2(czrg,:,:,:);
            % calculate the QUV from the two-channel complex
            for iR = 1:nr
                for iL = 1: nWin
                    if gpu_available
                        [S0(:,:,iR,iL),S1(:,:,iR,iL),S2(:,:,iR,iL),S3(:,:,iR,iL)] = ...
                            cumulativeQUV_gpu(IMGs_ch1(:,:,iR,iL),IMGs_ch2(:,:,iR,iL));
                    else
                        % 内联实现cumulativeQUV功能，避免调用嵌套函数
                        axis=angle(IMGs_ch2(:,:,iR,iL).*conj(IMGs_ch1(:,:,iR,iL)));
                        S0(:,:,iR,iL)=abs(IMGs_ch1(:,:,iR,iL)).^2+abs(IMGs_ch2(:,:,iR,iL)).^2;
                        S1(:,:,iR,iL)=abs(IMGs_ch1(:,:,iR,iL)).^2-abs(IMGs_ch2(:,:,iR,iL)).^2;
                        S2(:,:,iR,iL)=2.*abs(IMGs_ch1(:,:,iR,iL)).*abs(IMGs_ch2(:,:,iR,iL)).*cos(axis);
                        S3(:,:,iR,iL)=2.*abs(IMGs_ch1(:,:,iR,iL)).*abs(IMGs_ch2(:,:,iR,iL)).*sin(-axis);
                    end
                end
            end
            % average QUV accross split-spectrum and calculate dopu 
            dopu_ss = sqrt(mean(S1./S0,4).^2+mean(S2./S0,4).^2+mean(S3./S0,4).^2);%% ===>nZcrop,nX,nr
            dopu_splitSpectrum(:,:,iY) = gather(mean(dopu_ss,3));
        end
        
        %% whole spectrum nZ*nX*nR==> fft(complex)
        Bimg1_wholeStr=fft(Bd1.*gpu_winG_whole,SPL,1);Bimg2_wholeStr=fft(Bd2.*gpu_winG_whole,SPL,1);
        IMG1_wholeStr = Bimg1_wholeStr(czrg,:,:);IMG2_wholeStr = Bimg2_wholeStr(czrg,:,:);

        %% Struc Stokes, and OAC
        % 直接在循环内实现cumulativeQUV功能，避免调用嵌套函数
        if gpu_available
            [wS0,wS1,wS2,wS3] = cumulativeQUV_gpu(IMG1_wholeStr,IMG2_wholeStr);
        else
            % 内联实现cumulativeQUV功能
            axis=angle(IMG2_wholeStr.*conj(IMG1_wholeStr));
            wS0=abs(IMG1_wholeStr).^2+abs(IMG2_wholeStr).^2;
            wS1=abs(IMG1_wholeStr).^2-abs(IMG2_wholeStr).^2;
            wS2=2.*abs(IMG1_wholeStr).*abs(IMG2_wholeStr).*cos(axis);
            wS3=2.*abs(IMG1_wholeStr).*abs(IMG2_wholeStr).*sin(-axis);
        end
        wQ = wS1./wS0;wU = wS2./wS0; wV = wS3./wS0;
        strLin = mean(wS0,3);
        Strus(:,:,iY) = gather(20*log10(strLin));
        Smap_avg(:,:,:,iY) = gather(cat(3,mean(wQ,3),mean(wU,3),mean(wV,3)));
        Smap_rep1(:,:,:,iY) = gather(cat(3,wQ(:,:,1),wU(:,:,1),wV(:,:,1)));
        
        if ~hasSeg
            % 直接在循环内实现calOAC功能，避免调用嵌套函数
            % 将strLin转移到CPU进行处理
            cpu_strLin = gather(strLin);
            local_OAC = cpu_strLin*0;
            Nois_M = mean(mean(cpu_strLin(size(cpu_strLin,1)-2:size(cpu_strLin,1),:)));
            if Nois_M >= 1
                Nois_D = std(std(cpu_strLin(size(cpu_strLin,1)-2:size(cpu_strLin,1),:)));
                cpu_strLin = cpu_strLin-Nois_M+Nois_D; %去除噪声地板
                down = cpu_strLin<0;
                cpu_strLin(down)=1; %去除负数部分
                tail_signal = mean(cpu_strLin(size(cpu_strLin,1)-5:size(cpu_strLin,1),:)); %估计剩余光强
                tail_signal = medfilt1(tail_signal,15);
                tail_signal = smooth(tail_signal, 15)';

                for z=1:size(cpu_strLin,1)
                    local_OAC(z,:) = cpu_strLin(z,:)./(2*0.0086*sum(cpu_strLin(z+1:size(cpu_strLin,1),:))+tail_signal);
                end
            end
            % 使用外部函数surf_seg
            if exist('surf_seg.m', 'file') == 2
                topLines(:,iY) = surf_seg(local_OAC,0.25)+2;
            else
                fprintf('警告: surf_seg.m 函数未找到，无法计算表面分割\n');
            end
            Stru_OAC(:,:,iY) = local_OAC;
        end
        
        %% drLA and drPhR
        dopu_splitSpec_M = squeeze(mean(dopu_ss,3)); %% dopu across the different repeat(nR)
        if do_cfg1 
            if do_avg %% avg -- cfg1
                IMG_ch1 = squeeze(mean(IMG1_wholeStr,3));IMG_ch2 = squeeze(mean(IMG2_wholeStr,3));
                % 检查函数是否存在
                if exist('calLAPhRALL_gpu.m', 'file') == 2
                    [LA_c_cfg1_avg(:,:,:,iY),PhR_c_cfg1_avg(:,:,iY),cumLA_cfg1_avg(:,:,:,iY),...
                        LA_Ms_cfg1_rmBG(:,:,:,iY),PhR_Ms_cfg1_rmBG(:,:,iY),cumLA_Ms_cfg1_rmBG(:,:,:,iY)] = ...
                        calLAPhRALL_gpu(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M,kRL_cfg1,kRU_cfg1,h1,h2,Avnum,wovWinF,gpu_available);
                else
                    fprintf('警告: calLAPhRALL_gpu.m 函数未找到，跳过 cfg1 avg 计算\n');
                end
            end
            if do_eig %% eig -- cfg1
                % 检查函数是否存在
                if exist('OCTA_F_ED_Clutter_EigFeed_comp.m', 'file') == 2 && exist('calLAPhRALL_gpu.m', 'file') == 2
                    [IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,IMG1_wholeStr,IMG2_wholeStr), 1); % 提取静态分量
                    IMG_ch1 = IMG_ch(1:nZ,:);IMG_ch2 = IMG_ch(nZ+1:end,:);
                    [LA_c_cfg1_eig(:,:,:,iY),PhR_c_cfg1_eig(:,:,iY),cumLA_cfg1_eig(:,:,:,iY)] = ...
                    calLAPhRALL_gpu(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M,kRL_cfg1,kRU_cfg1,h1,h2,Avnum,wovWinF,gpu_available);
                else
                    if exist('OCTA_F_ED_Clutter_EigFeed_comp.m', 'file') ~= 2
                        fprintf('警告: OCTA_F_ED_Clutter_EigFeed_comp.m 函数未找到，跳过 cfg1 eig 计算\n');
                    end
                    if exist('calLAPhRALL_gpu.m', 'file') ~= 2
                        fprintf('警告: calLAPhRALL_gpu.m 函数未找到，跳过 cfg1 eig 计算\n');
                    end
                end
            end
        end
        if do_cfg2 
            if do_avg %% avg -- cfg2
                IMG_ch1 = squeeze(mean(IMG1_wholeStr,3));IMG_ch2 = squeeze(mean(IMG2_wholeStr,3));
                % 检查函数是否存在
                if exist('calLAPhRcfg2_gpu.m', 'file') == 2
                    [LA_c_cfg2_avg(:,:,:,iY),PhR_c_cfg2_avg(:,:,iY),cumLA_cfg2_avg(:,:,:,iY)] = ...
                        calLAPhRcfg2_gpu(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M,kRL_cfg2,kRU_cfg2,h1,h2,Avnum,wovWinF,gpu_available);
                else
                    fprintf('警告: calLAPhRcfg2_gpu.m 函数未找到，跳过 cfg2 avg 计算\n');
                end
            end
            if do_eig %% eig -- cfg2
                % 检查函数是否存在
                if exist('OCTA_F_ED_Clutter_EigFeed_comp.m', 'file') == 2 && exist('calLAPhRcfg2_gpu.m', 'file') == 2
                    [IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,IMG1_wholeStr,IMG2_wholeStr), 1); % 提取静态分量
                    IMG_ch1 = IMG_ch(1:nZ,:);IMG_ch2 = IMG_ch(nZ+1:end,:);
                    [LA_c_cfg2_eig(:,:,:,iY),PhR_c_cfg2_eig(:,:,iY),cumLA_cfg2_eig(:,:,:,iY)] = ...
                    calLAPhRcfg2_gpu(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M,kRL_cfg2,kRU_cfg2,h1,h2,Avnum,wovWinF,gpu_available);
                else
                    if exist('OCTA_F_ED_Clutter_EigFeed_comp.m', 'file') ~= 2
                        fprintf('警告: OCTA_F_ED_Clutter_EigFeed_comp.m 函数未找到，跳过 cfg2 eig 计算\n');
                    end
                    if exist('calLAPhRcfg2_gpu.m', 'file') ~= 2
                        fprintf('警告: calLAPhRcfg2_gpu.m 函数未找到，跳过 cfg2 eig 计算\n');
                    end
                end
            end
        end 
        %% check results from a frame
        if show_img
            imgcRg = 1:320;
            figure(1);imshow(20*log10(abs(IMG2_wholeStr(:,:,1))),[80 120]);title('single channel rep1 Structure')
            figure(2);imshow(20*log10(mean(abs(IMG1_wholeStr(:,:,:)),3)),[80 120]);title('single channel avg Structure')
            figure(3);imshow(Strus(:,:,iY),[170 245]);title('merged avg Stru')
            figure(4);imshow(Smap_rep1(:,:,:,iY),[]);title('rep1 Stokes')
            figure(5);imshow(Smap_avg(:,:,:,iY),[]);title('avg Stokes')
            figure(6);imshow(local_OAC,[0 3]);title('oac')
            hold on;plot(topLines(:,iY),'-r');hold off;
            
            if do_cfg2 && do_avg
                figure(7);imshow(LA_c_cfg2_avg(:,:,:,iY)/2+0.5,[]);title('LA cfg2 avg')
                figure(8);imagesc(PhR_c_cfg2_avg(:,:,iY),[0 0.5]);title('PhR cfg2 avg')
                figure(9);imshow(cumLA_cfg2_avg(:,:,:,iY)/2+0.5,[]);title('cumLA cfg2 avg')
            end
            if do_cfg1 && do_eig
                figure(10);imshow(LA_c_cfg1_eig(:,:,:,iY)/2+0.5,[]);title('LA cfg1 eig')
                figure(11);imagesc(PhR_c_cfg1_eig(:,:,iY),[0 0.5]);title('PhR cfg1 eig')
            end
        end
        fclose(fid);
        pb.progress;
    end
    pb.stop;
    fclose all;
    
    % 从cell数组中合并结果到全局数组
    % 已经在for循环中直接存储到全局数组，不需要合并
    
    %% save results: strus(flow),stokes,oac
    if saveDicom
        % 修复索引越界问题
        slice_index = min(100, floor(size(Strus,3)/2));  % 如果不足100层，则取中间层
        if slice_index < 1
            slice_index = 1;  % 至少取第一层
        end
        SS=Strus(:,:,slice_index);strUrg = max(SS(:))-5;strLrg = min(SS(:))+5;
        
        % 预分配Struc数组
        Struc = zeros(size(Strus,1), size(Strus,2), 1, size(Strus,3));
        
        for i=1:size(Strus,3)
            SS1=Strus(:,:,i);
            Struc(:,:,:,i)=(SS1-strLrg)./(strUrg-strLrg);
        end
        
        [Struc_flat] = volFlatten(Struc,topLines);
        dicomwrite(uint8(255*(Struc)),[foutputdir,name,'_1-1_Struc.dcm']);
        dicomwrite(uint8(255*(Struc_flat)),[foutputdir,name,'_1-1_Struc_flat.dcm']);
        
        writematrix([strLrg strUrg],[foutputdir,name,'_1-1_StrucRg.txt']);
        dicomwrite(uint8(255*(Smap_rep1/2+0.5)),[foutputdir,name,'_1-3_1rep-Stokes.dcm']);
        dicomwrite(uint8(255*(Smap_avg/2+0.5)),[foutputdir,name,'_1-3_4rep-Stokes.dcm']);
        dicomwrite(uint8(255*(permute(dopu_splitSpectrum,[1 2 4 3]))),[foutputdir,name,'_1-4_dopu_SS.dcm']);

        if ~hasSeg
            save([foutputdir,name, '.mat'],'topLines','czrg');
        end
        save([foutputdir,'topLines.mat'],'topLines','czrg');
        rotAngle = 440;
        
        if do_cfg1
            PRRrg = [0 0.5];
            writematrix(PRRrg,[foutputdir,name,'_2-0_PhRRg.txt']);
            if do_avg
                % 预分配变量
                PRRc = zeros(size(PhR_c_cfg1_avg,1), size(PhR_c_cfg1_avg,2), 3, nY, 'uint8');
                cumLA_cfg_hsv = zeros(size(cumLA_cfg1_avg,1), size(cumLA_cfg1_avg,2), 3, nY);
                LA_cfg_hsv = zeros(size(LA_c_cfg1_avg,1), size(LA_c_cfg1_avg,2), 3, nY);
                PRRc_rmBG = zeros(size(PhR_Ms_cfg1_rmBG,1), size(PhR_Ms_cfg1_rmBG,2), 3, nY, 'uint8');
                cumLA_Ms_cfg1_rmBG_hsv = zeros(size(cumLA_Ms_cfg1_rmBG,1), size(cumLA_Ms_cfg1_rmBG,2), 3, nY);
                LA_Ms_cfg1_rmBG_hsv = zeros(size(LA_Ms_cfg1_rmBG,1), size(LA_Ms_cfg1_rmBG,2), 3, nY);
                
                for iY =1:nY
                    PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_avg(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg1_avg(:,:,:,iY),rotAngle);
                    [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg1_avg(:,:,:,iY),rotAngle);
                    
                    PRRc_rmBG(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_Ms_cfg1_rmBG(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_Ms_cfg1_rmBG_hsv(:,:,:,iY)] = quColoring(cumLA_Ms_cfg1_rmBG(:,:,:,iY),rotAngle);
                    [LA_Ms_cfg1_rmBG_hsv(:,:,:,iY)] = quColoring(LA_Ms_cfg1_rmBG(:,:,:,iY),rotAngle);
                end
                
                dicomwrite(uint8(255*(cumLA_cfg1_avg/2+0.5)),[foutputdir,name,'_2-1_cumLA-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_2-2_cumLA-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_c_cfg1_avg/2+0.5)),[foutputdir,name,'_2-3_drLA-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_2-4_drLA-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(PRRc,[foutputdir,name,'_2-5_PhR-',num2str(nr),'repAvg.dcm']);
                
                dicomwrite(uint8(255*(cumLA_Ms_cfg1_rmBG/2+0.5)),[foutputdir,name,'_2-6_cumLA_rmBG-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*cumLA_Ms_cfg1_rmBG_hsv),[foutputdir,name,'_2-7_cumLA_rmBG-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_Ms_cfg1_rmBG/2+0.5)),[foutputdir,name,'_2-8_drLA_rmBG-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*LA_Ms_cfg1_rmBG_hsv),[foutputdir,name,'_2-9_drLA_rmBG-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(PRRc_rmBG,[foutputdir,name,'_2-10_PhR_rmBG-',num2str(nr),'repAvg.dcm']);
            end
            
            if do_eig
                % 预分配变量
                PRRc = zeros(size(PhR_c_cfg1_eig,1), size(PhR_c_cfg1_eig,2), 3, nY, 'uint8');
                cumLA_cfg_hsv = zeros(size(cumLA_cfg1_eig,1), size(cumLA_cfg1_eig,2), 3, nY);
                LA_cfg_hsv = zeros(size(LA_c_cfg1_eig,1), size(LA_c_cfg1_eig,2), 3, nY);
                
                for iY =1:nY
                    PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_eig(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg1_eig(:,:,:,iY),rotAngle);
                    [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg1_eig(:,:,:,iY),rotAngle);
                end
                
                dicomwrite(uint8(255*(cumLA_cfg1_eig/2+0.5)),[foutputdir,name,'_2-6_cumLA-',num2str(nr),'repEig.dcm']);
                dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_2-7_cumLA-',num2str(nr),'repEig_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_c_cfg1_eig/2+0.5)),[foutputdir,name,'_2-8_drLA-',num2str(nr),'repEig.dcm']);
                dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_2-9_drLA-',num2str(nr),'repEig_hsvColoring.dcm']);
                dicomwrite(PRRc,[foutputdir,name,'_2-10_PhR-',num2str(nr),'repEig.dcm']);
            end
        end

        if do_cfg2
            PRRrg = [0 0.5];
            writematrix(PRRrg,[foutputdir,name,'_3-0_PhRRg.txt']);
            
            if do_avg
                % 预分配变量
                PRRc = zeros(size(PhR_c_cfg2_avg,1), size(PhR_c_cfg2_avg,2), 3, nY, 'uint8');
                cumLA_cfg_hsv = zeros(size(cumLA_cfg2_avg,1), size(cumLA_cfg2_avg,2), 3, nY);
                LA_cfg_hsv = zeros(size(LA_c_cfg2_avg,1), size(LA_c_cfg2_avg,2), 3, nY);
                
                for iY =1:nY
                    PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg2_avg(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg2_avg(:,:,:,iY));
                    [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg2_avg(:,:,:,iY));
                end
                
                dicomwrite(uint8(255*(cumLA_cfg2_avg/2+0.5)),[foutputdir,name,'_3-1_cumLA-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_3-2_cumLA-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_c_cfg2_avg/2+0.5)),[foutputdir,name,'_3-3_drLA-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_3-4_drLA-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(PRRc,[foutputdir,name,'_3-5_PhR-',num2str(nr),'repAvg.dcm']);
            end
            
            if do_eig
                % 预分配变量
                PRRc = zeros(size(PhR_c_cfg2_eig,1), size(PhR_c_cfg2_eig,2), 3, nY, 'uint8');
                cumLA_cfg_hsv = zeros(size(cumLA_cfg2_eig,1), size(cumLA_cfg2_eig,2), 3, nY);
                LA_cfg_hsv = zeros(size(LA_c_cfg2_eig,1), size(LA_c_cfg2_eig,2), 3, nY);
                
                for iY =1:nY
                    PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg2_eig(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg2_eig(:,:,:,iY));
                    [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg2_eig(:,:,:,iY));
                end
                
                dicomwrite(uint8(255*(cumLA_cfg2_eig/2+0.5)),[foutputdir,name,'_3-6_cumLA-',num2str(nr),'repEig.dcm']);
                dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_3-7_cumLA-',num2str(nr),'repEig_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_c_cfg2_eig/2+0.5)),[foutputdir,name,'_3-8_drLA-',num2str(nr),'repEig.dcm']);
                dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_3-9_drLA-',num2str(nr),'repEig_hsvColoring.dcm']);
                dicomwrite(PRRc,[foutputdir,name,'_3-10_PhR-',num2str(nr),'repEig.dcm']);
            end
        end
    end % save results
    
    fprintf('文件 %s 处理完成!\n', display_name);
    fprintf('输出目录: %s\n', foutputdir);
    fprintf('\n');
end % for nr = [nR]

%% 释放GPU内存
if gpu_available
    reset(gpuDevice);
end

end % function rPSOCT_process_single_file_gpu