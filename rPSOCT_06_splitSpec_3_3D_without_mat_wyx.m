%%% rPSOCT_process_single_file.m
%%% 1. frame based code to cal 3D LAs and Phase[two cfg: cfg1== LA, cfg2==Phase]
%%% 1. use doup to apply a variable gaussian fitler to Q U V before cal LA
%%%   and phase
%%% 2. speed the cal process, setZrg = 100; setZrg = 0 for process all
%%% ______________________________20240913________________________________________                if do_avg %% avg -- cfg1
                    if gpu_available
                        IMG_ch1 = gather(squeeze(mean(IMG1_wholeStr,3)));
                        IMG_ch2 = gather(squeeze(mean(IMG2_wholeStr,3)));
                        % 处理dopu_ss的不同情况
                        if ~do_ssdopu
                            dopu_splitSpec_M_cpu = 1;
                        else
                            dopu_splitSpec_M_cpu = gather(squeeze(mean(dopu_ss,3)));
                        end
                    else
                        IMG_ch1 = squeeze(mean(IMG1_wholeStr,3));
                        IMG_ch2 = squeeze(mean(IMG2_wholeStr,3));
                        % 处理dopu_ss的不同情况
                        if ~do_ssdopu
                            dopu_splitSpec_M_cpu = 1;
                        else
                            dopu_splitSpec_M_cpu = squeeze(mean(dopu_ss,3));
                        end
                    end%% 3. optional to toggle vWinF
%%% 4.cfg2 was optinal to [do SVD with remove bias] and DDG 
%%% ______________________________20240721_______________________________________________
%%% 5. optinal to perform split-spectrum DOPU
%%% 6. DIR: 03_2024.07.03.psOCTry2

% 添加function文件夹到搜索路径
script_dir = fileparts(mfilename('fullpath'));
function_path = fullfile(script_dir, 'function');
if exist(function_path, 'dir')
    addpath(function_path);
end

% 设置数据路径
data_path = 'D:\1-Liu Jian\yongxin.wang\PSOCT\data_without_mat\';
output_base = 'D:\1-Liu Jian\yongxin.wang\PSOCT\without_wovWinf\';

% 检查输入路径
if ~exist(data_path, 'dir')
    error(['数据路径不存在: ' data_path]);
end

% 创建输出路径
if ~exist(output_base, 'dir')
    mkdir(output_base);
end

% 获取所有oct文件
disp('正在处理数据...');
oct_files = dir(fullfile(data_path, '*.oct'));

if isempty(oct_files)
else
% 显示找到的文件
fprintf('找到 %d 个 OCT 文件:\n', length(oct_files));
for i = 1:length(oct_files)
    fprintf('[%d] %s\n', i, oct_files(i).name);
end

% 逐个处理文件
for i = 1:length(oct_files)
    full_path = fullfile(data_path, oct_files(i).name);
    fprintf('\n==================================================\n');
    fprintf('开始处理第 %d/%d 个文件: %s\n', i, length(oct_files), oct_files(i).name);
    fprintf('==================================================\n');
    
    try
        % 调用单文件处理函数
        rPSOCT_process_single_file(full_path, output_base);
        fprintf('文件 %s 处理成功\n', oct_files(i).name);
    catch ME
        fprintf('处理文件 %s 时出错: %s\n', oct_files(i).name, ME.message);
        % 继续处理下一个文件
    end
end

fprintf('\n所有文件处理完成!\n');

end


function rPSOCT_process_single_file(varargin)
    %% (1) 获取输入参数
    if nargin < 1
        error('请提供输入文件路径作为参数');
    end
    
    input_file_path = varargin{1};
    
    % 如果提供了输出路径，则使用提供的路径，否则使用默认路径
    if nargin >= 2
        output_base = varargin{2};
    else
        output_base = 'D:\1-Liu Jian\yongxin.wang\DOPU_output\';
    end
    
    % 检查输入文件是否存在
    if ~exist(input_file_path, 'file')
        error(['输入文件不存在: ' input_file_path]);
    end
    
    % 创建输出目录
    if ~exist(output_base,'dir'), mkdir(output_base), end
    
    %% (2) 初始化参数
    %clear all;close all;clc;
    delete(gcp('nocreate'));
    
    % 检测GPU可用性
    try
        gpu_info = gpuDevice();
        gpu_available = true;
        fprintf('检测到GPU: %s (内存: %.2f GB)\n', gpu_info.Name, gpu_info.AvailableMemory/1024^3);
        fprintf('GPU计算已启用\n');
    catch
        gpu_available = false;
        fprintf('未检测到支持的GPU或CUDA，将使用CPU计算\n');
    end
    
    % 启动并行池
    if isempty(gcp('nocreate'))
        parpool('local', 34);
    end
    
    disp_coef=-20.1;
    do_PhComp=1; do_medianshift=1;
    do_reg=0; % toggle Inter-Frame registration (for motion removal)
    useref=1;
    saveDicom=1;
    show_img=0;
    iy=1;
    hasSeg = 0;%is have mat?

    distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing

    %% (3) 处理参数设置
    do_cfg1 = 1;
    do_cfg2 = 0;
    do_avg = 1;
    do_eig =0;
    setZrg = 0; %% if speed the cal process setZrg = 200 to speed
    wovWinF = 1; % if do variabe weight filter use DOPU

    %% 设置是否采用ssdopu 
    do_ssdopu = 0; % if do splitted spectrum
    %% (3)偏振参数设置 drLA and PhR
    Avnum = 7;   %number of layers made use in first calculation to compute the normal vector in Poincare sphere
    kRL_cfg1 = 2;kRU_cfg1 = 7;% for LA calculation
    kRL_cfg2 = 3;kRU_cfg2 = 21;% for phase retardation calculation
    Avnum = 3; % for DDG test
    h1 = fspecial('gaussian',[21 21],9);%预处理 corr3 
    h2 = fspecial('gaussian',[17 17],9);
    
    %% （4）处理单个文件
    filename = input_file_path;
    [filepath, name, ext] = fileparts(filename);
    display_name = [name ext];
    
    fprintf('\n==================================================\n');
    fprintf('正在处理文件: %s\n', display_name);
    fprintf('文件路径: %s\n', filename);
    fprintf('==================================================\n');
    
    Thr=170; % Threshold to remove low OCT signal
    % This reads the parameter used for the data acquisition from *.oct* file
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
    K=1:Blength;
    use_autoRg=1;
    RgFlow=[60 110];
    jsurf=zeros(4);
    IMcropRg=1:(SPL/2);
    nZcrop=numel(IMcropRg);
    imshowrgZ=1:nZcrop;
    jusamp=zeros(nZcrop,nX,4);
    winG =   tukeywin(Blength,0.25);
    %(a)subWins: params for split spectrum DOPU
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

   
%(b) Calculates the Cumulative Stokes parameters I,Q,U,V 
for nr = [nR]
    if setZrg
        foutputdir = fullfile(output_base,[name,'.rep',num2str(nr),'fastZ_testDDG',num2str(setZrg),'\']);
    else
        foutputdir = fullfile(output_base,[name,'.rep',num2str(nr),'DDG_corr3_wobg_test\']);
    end
    if ~exist(foutputdir,'dir'), mkdir(foutputdir); end
    
    czrg = [21:320];% set z range
    topLines =ones(nX,nY);
    if hasSeg
        load([name, '.mat']);
        topLines = round(topLines);topLines(topLines<=1) = 1;
    end
    if setZrg, czrg = czrg(1:setZrg);end
    nZcrop = numel(czrg);nZ = nZcrop;
    % (c) creat array to store results
    LA_c_cfg1_avg = zeros(nZcrop-Avnum,nX,3,nY);PhR_c_cfg1_avg = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg1_avg = zeros(nZcrop-Avnum,nX,3,nY);
    LA_c_cfg1_eig = zeros(nZcrop-Avnum,nX,3,nY);PhR_c_cfg1_eig = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg1_eig = zeros(nZcrop-Avnum,nX,3,nY);

    LA_Ms_cfg1_rmBG = LA_c_cfg1_avg;
    PhR_Ms_cfg1_rmBG = cumLA_cfg1_avg;
    cumLA_Ms_cfg1_rmBG = LA_c_cfg1_eig;

    LA_c_cfg2_avg = zeros(nZcrop-Avnum,nX,3,nY);PhR_c_cfg2_avg = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg2_avg = zeros(nZcrop-Avnum,nX,3,nY);
    LA_c_cfg2_eig = zeros(nZcrop-Avnum,nX,3,nY);PhR_c_cfg2_eig = zeros(nZcrop-Avnum,nX,nY);
    cumLA_cfg2_eig = zeros(nZcrop-Avnum,nX,3,nY);
    Smap_avg = zeros(numel(czrg),nX,3,nY);Smap_rep1 = zeros(numel(czrg),nX,3,nY);
    Strus = zeros(numel(czrg),nX,nY);Stru_OAC = Strus;
    dopu_splitSpectrum = zeros(numel(czrg),nX,nY);
    pb = ProgressBar(nY);

    parfor iY=1:nY-1 % Bscan number
    % for iY=1:50:nY-1 % Bscan number
    % for iY = 478 % ===============================>>>>>>>>>>> run single B scan
            
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
                Bd1 =real(hilbert(Bs1).*phV);  % dispersion correction
                Bd2 =real(hilbert(Bs2).*phV);
            else
                Bd1=Bs1;
                Bd2=Bs2;
            end
            
            % 如果GPU可用，将数据转移到GPU进行FFT计算
            if gpu_available
                Bd1_gpu = gpuArray(Bd1);
                Bd2_gpu = gpuArray(Bd2);
                winG_gpu = gpuArray(winG);
                winG_whole_gpu = gpuArray(winG_whole);
            else
                Bd1_gpu = Bd1;
                Bd2_gpu = Bd2;
                winG_gpu = winG;
                winG_whole_gpu = winG_whole;
            end
            %% split spectrum
            if ~do_ssdopu % if not do adaptive gaussian filter dopu_ss =1 
                dopu_ss = 1;
                % 设置默认值，但在后续处理中需要正确设置
            else
                % creat array to store split-spectrum complex (FFT): Z*X*nR*nWin 
                if gpu_available
                    Bimg1 = zeros(SPL,nX,nr,nWin,'gpuArray');Bimg2 = Bimg1;
                    S0 = zeros(nZcrop,nX,nr,nWin,'gpuArray'); S1 = S0;S2=S0;S3=S0;
                else
                    Bimg1 = zeros(SPL,nX,nr,nWin);Bimg2 = Bimg1;
                    S0 = zeros(nZcrop,nX,nr,nWin); S1 = S0;S2=S0;S3=S0;
                end
                
                for iR = 1:nr
                    for iL=1:nWin
                        % extract data fragments from two different channel and
                        % apply a giassian filter before performing FFT
                        iBd1=Bd1_gpu(windex(iL):windex(iL)+winL-1,:,iR).*winG_gpu;
                        iBd2=Bd2_gpu(windex(iL):windex(iL)+winL-1,:,iR).*winG_gpu;
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
                            [S0(:,:,iR,iL),S1(:,:,iR,iL),S2(:,:,iR,iL),S3(:,:,iR,iL)] = ...
                                cumulativeQUV(IMGs_ch1(:,:,iR,iL),IMGs_ch2(:,:,iR,iL));
                        end
                    end
                end
                % average QUV accross split-spectrum and calculate dopu 
                dopu_ss = sqrt(mean(S1./S0,4).^2+mean(S2./S0,4).^2+mean(S3./S0,4).^2);%% ===>nZcrop,nX,nr
                if gpu_available
                    dopu_splitSpectrum(:,:,iY) = gather(mean(dopu_ss,3));
                else
                    dopu_splitSpectrum(:,:,iY) = mean(dopu_ss,3);
                end
            end
            %% whole spectrum nZ*nX*nR==> fft(complex)
            Bimg1_wholeStr=fft(Bd1_gpu.*winG_whole_gpu,SPL,1);Bimg2_wholeStr=fft(Bd2_gpu.*winG_whole_gpu,SPL,1);
            IMG1_wholeStr = Bimg1_wholeStr(czrg,:,:);IMG2_wholeStr = Bimg2_wholeStr(czrg,:,:);

            %% Struc Stokes, and OAC
            if gpu_available
                [wS0,wS1,wS2,wS3] = cumulativeQUV_gpu(IMG1_wholeStr,IMG2_wholeStr);
            else
                [wS0,wS1,wS2,wS3] = cumulativeQUV(IMG1_wholeStr,IMG2_wholeStr);
            end
            wQ = wS1./wS0;wU = wS2./wS0; wV = wS3./wS0;
            strLin = mean(wS0,3);
            
            % 将GPU计算结果转回CPU以保存
            if gpu_available
                Strus(:,:,iY) = gather(20*log10(strLin));
                Smap_avg(:,:,:,iY) = gather(cat(3,mean(wQ,3),mean(wU,3),mean(wV,3)));
                Smap_rep1(:,:,:,iY) = gather(cat(3,wQ(:,:,1),wU(:,:,1),wV(:,:,1)));
                strLin_cpu = gather(strLin); % 用于后续OAC计算
            else
                Strus(:,:,iY) = 20*log10(strLin);
                Smap_avg(:,:,:,iY) = cat(3,mean(wQ,3),mean(wU,3),mean(wV,3));
                Smap_rep1(:,:,:,iY) = cat(3,wQ(:,:,1),wU(:,:,1),wV(:,:,1));
                strLin_cpu = strLin;
            end
            if ~hasSeg, strOAC = calOAC(strLin_cpu); topLines(:,iY)=surf_seg(strOAC,0.25)+2; end
            
            %% drLA and drPhR
            % dopu across the different repeat(nR) - now handled per GPU/CPU branch above
            if do_cfg1 
                if do_avg %% avg -- cfg1
                    if gpu_available
                        IMG_ch1 = gather(squeeze(mean(IMG1_wholeStr,3)));
                        IMG_ch2 = gather(squeeze(mean(IMG2_wholeStr,3)));
                        dopu_splitSpec_M_cpu = gather(squeeze(mean(dopu_ss,3)));
                    else
                        IMG_ch1 = squeeze(mean(IMG1_wholeStr,3));
                        IMG_ch2 = squeeze(mean(IMG2_wholeStr,3));
                        dopu_splitSpec_M_cpu = squeeze(mean(dopu_ss,3));
                    end
                    [LA_c_cfg1_avg(:,:,:,iY),PhR_c_cfg1_avg(:,:,iY),cumLA_cfg1_avg(:,:,:,iY),...
                        LA_Ms_cfg1_rmBG(:,:,:,iY),PhR_Ms_cfg1_rmBG(:,:,iY),cumLA_Ms_cfg1_rmBG(:,:,:,iY)] = ...
                        calLAPhRALL(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M_cpu,kRL_cfg1,kRU_cfg1,h1,h2,Avnum,wovWinF);
                end
                if do_eig %% eig -- cfg1
                    if gpu_available
                        [IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,gather(IMG1_wholeStr),gather(IMG2_wholeStr)), 1); % 提取静态分量
                    else
                        [IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,IMG1_wholeStr,IMG2_wholeStr), 1); % 提取静态分量
                    end
                    IMG_ch1 = IMG_ch(1:nZ,:);IMG_ch2 = IMG_ch(nZ+1:end,:);
                    [LA_c_cfg1_eig(:,:,:,iY),PhR_c_cfg1_eig(:,:,iY),cumLA_cfg1_eig(:,:,:,iY)] = ...
                    calLAPhRALL(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M_cpu,kRL_cfg1,kRU_cfg1,h1,h2,Avnum,wovWinF);
                end
            end
            if do_cfg2 
                if do_avg %% avg -- cfg2
                    if gpu_available
                        IMG_ch1 = gather(squeeze(mean(IMG1_wholeStr,3)));
                        IMG_ch2 = gather(squeeze(mean(IMG2_wholeStr,3)));
                    else
                        IMG_ch1 = squeeze(mean(IMG1_wholeStr,3));
                        IMG_ch2 = squeeze(mean(IMG2_wholeStr,3));
                    end
                    [LA_c_cfg2_avg(:,:,:,iY),PhR_c_cfg2_avg(:,:,iY),cumLA_cfg2_avg(:,:,:,iY)] = ...
                        calLAPhRcfg2(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M_cpu,kRL_cfg2,kRU_cfg2,h1,h2,Avnum,wovWinF);
                end
                if do_eig %% eig -- cfg2
                    if gpu_available
                        [IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,gather(IMG1_wholeStr),gather(IMG2_wholeStr)), 1); % 提取静态分量
                    else
                        [IMG_ch,~] = OCTA_F_ED_Clutter_EigFeed_comp(cat(1,IMG1_wholeStr,IMG2_wholeStr), 1); % 提取静态分量
                    end
                    IMG_ch1 = IMG_ch(1:nZ,:);IMG_ch2 = IMG_ch(nZ+1:end,:);
                    [LA_c_cfg2_eig(:,:,:,iY),PhR_c_cfg2_eig(:,:,iY),cumLA_cfg2_eig(:,:,:,iY)] = ...
                    calLAPhRcfg2(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M_cpu,kRL_cfg2,kRU_cfg2,h1,h2,Avnum,wovWinF);
                end
            end 
            %% check results from a frame
            if show_img
            imgcRg = 1:320;
    %         figure(1);imshow(20*log10(abs(IMG2_wholeStr(:,:,1))),[80 120]);title('single channel rep1 Structure')
    %         figure(2);imshow(20*log10(mean(abs(IMG1_wholeStr(:,:,:)),3)),[80 120]);title('single channel avg Structure')
    %         figure(3);imshow(squeeze(Strus(:,:,iY)),[170 245]);title('merged avg Stru')
    %         figure(4);imshow(squeeze(Smap_rep1(:,:,:,iY)),[]);title('rep1 Stokes')
    %         figure(5);imshow(squeeze(Smap_avg(:,:,:,iY)),[]);title('avg Stokes')
    % %         figure(6);imshow(squeeze(Stru_OAC(:,:,iY)),[0 3]);title('oac')
    %         hold on;plot(topLines(:,iY),'-r');hold off;
    %         figure(7);imshow(squeeze(LA_c_cfg2_avg(:,:,:,iY))/2+0.5,[]);title('LA cfg1 avg')
    %         figure(8);imagesc(squeeze(PhR_c_cfg2_avg(:,:,200)),[0 0.5]);title('PhR cfg2 avg')
    %         figure(7);imshow(squeeze(cumLA_cfg2_avg(:,:,:,200))/2+0.5,[]);title('LA cfg2 avg')
    %         figure(9);imshow(squeeze(LA_c_cfg1_eig(:,:,:,iY))/2+0.5,[]);title('LA cfg1 eig')
    %         figure(10);imagesc(squeeze(PhR_c_cfg1_eig(:,:,iY)),[0 0.5]);title('PhR cfg1 eig')
            end
            
    %%
            fclose(fid);
            pb.progress;
    end
    pb.stop;
    fclose all;
    
    % 释放GPU内存
    if gpu_available
        reset(gpuDevice);
        fprintf('GPU内存已释放\n');
    end
    %% save results: strus(flow),stokes,oac
    if saveDicom
        % 修复索引越界问题
        slice_index = min(100, floor(size(Strus,3)/2));  % 如果不足100层，则取中间层
        if slice_index < 1
            slice_index = 1;  % 至少取第一层
        end
        SS=Strus(:,:,slice_index);strUrg = max(SS(:))-5;strLrg = min(SS(:))+5;
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
            for iY =1:nY
                PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_avg(:,:,iY),PRRrg)*256),parula(256))*256);
                [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg1_avg(:,:,:,iY),rotAngle);
                [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg1_avg(:,:,:,iY),rotAngle);
                
                PRRc_rmBG(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_Ms_cfg1_rmBG(:,:,iY),PRRrg)*256),parula(256))*256);
                [cumLA_Ms_cfg1_rmBG_hsv(:,:,:,iY)] = quColoring(cumLA_Ms_cfg1_rmBG(:,:,:,iY),rotAngle);
                [LA_Ms_cfg1_rmBG_hsv(:,:,:,iY)] = quColoring(LA_Ms_cfg1_rmBG(:,:,:,iY),rotAngle);
                
            end
            dicomwrite(uint8(255*(cumLA_cfg1_avg/2+0.5)),[foutputdir,name,'_2-1_cumLA-cfg1-',num2str(nr),'repAvg.dcm']);
            dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_2-2_cumLA-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
            dicomwrite(uint8(255*(LA_c_cfg1_avg/2+0.5)),[foutputdir,name,'_2-3_drLA-cfg1-',num2str(nr),'repAvg.dcm']);
            dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_2-4_drLA-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
            dicomwrite(PRRc,[foutputdir,name,'_2-5_PhR-cfg1-',num2str(nr),'repAvg.dcm']);
            
            dicomwrite(uint8(255*(cumLA_Ms_cfg1_rmBG/2+0.5)),[foutputdir,name,'_2-6_cumLA_rmBG-cfg1-',num2str(nr),'repAvg.dcm']);
            dicomwrite(uint8(255*cumLA_Ms_cfg1_rmBG_hsv),[foutputdir,name,'_2-7_cumLA_rmBG-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
            dicomwrite(uint8(255*(LA_Ms_cfg1_rmBG/2+0.5)),[foutputdir,name,'_2-8_drLA_rmBG-cfg1-',num2str(nr),'repAvg.dcm']);
            dicomwrite(uint8(255*LA_Ms_cfg1_rmBG_hsv),[foutputdir,name,'_2-9_drLA_rmBG-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
            dicomwrite(PRRc_rmBG,[foutputdir,name,'_2-10_PhR_rmBG-cfg1-',num2str(nr),'repAvg.dcm']);
        end
        if do_eig
            for iY =1:nY
                PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_eig(:,:,iY),PRRrg)*256),parula(256))*256);
                [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg1_eig(:,:,:,iY),rotAngle);
                [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg1_eig(:,:,:,iY),rotAngle);
            end
            dicomwrite(uint8(255*(cumLA_cfg1_eig/2+0.5)),[foutputdir,name,'_2-6_cumLA-cfg1-',num2str(nr),'repEig.dcm']);
            dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_2-7_cumLA-cfg1-',num2str(nr),'repEig_hsvColoring.dcm']);
            dicomwrite(uint8(255*(LA_c_cfg1_eig/2+0.5)),[foutputdir,name,'_2-8_drLA-cfg1-',num2str(nr),'repEig.dcm']);
            dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_2-9_drLA-cfg1-',num2str(nr),'repEig_hsvColoring.dcm']);
            dicomwrite(PRRc,[foutputdir,name,'_2-10_PhR-cfg1-',num2str(nr),'repEig.dcm']);
        end
        end

        if do_cfg2
        PRRrg = [0 0.5];
        writematrix(PRRrg,[foutputdir,name,'_3-0_PhRRg.txt']);
        if do_avg
            for iY =1:nY
                PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg2_avg(:,:,iY),PRRrg)*256),parula(256))*256);
                [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg2_avg(:,:,:,iY));
                [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg2_avg(:,:,:,iY));
            end
            dicomwrite(uint8(255*(cumLA_cfg2_avg/2+0.5)),[foutputdir,name,'_3-1_cumLA-cfg2-',num2str(nr),'repAvg.dcm']);
            dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_3-2_cumLA-cfg2-',num2str(nr),'repAvg_hsvColoring.dcm']);
            dicomwrite(uint8(255*(LA_c_cfg2_avg/2+0.5)),[foutputdir,name,'_3-3_drLA-cfg2-',num2str(nr),'repAvg.dcm']);
            dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_3-4_drLA-cfg2-',num2str(nr),'repAvg_hsvColoring.dcm']);
            dicomwrite(PRRc,[foutputdir,name,'_3-5_PhR-cfg2-',num2str(nr),'repAvg.dcm']);
        end
        if do_eig
            for iY =1:nY
                PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg2_eig(:,:,iY),PRRrg)*256),parula(256))*256);
                [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg2_eig(:,:,:,iY));
                [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg2_eig(:,:,:,iY));
            end
            dicomwrite(uint8(255*(cumLA_cfg2_eig/2+0.5)),[foutputdir,name,'_3-6_cumLA-cfg2-',num2str(nr),'repEig.dcm']);
            dicomwrite(uint8(255*cumLA_cfg_hsv),[foutputdir,name,'_3-7_cumLA-cfg2-',num2str(nr),'repEig_hsvColoring.dcm']);
            dicomwrite(uint8(255*(LA_c_cfg2_eig/2+0.5)),[foutputdir,name,'_3-8_drLA-cfg2-',num2str(nr),'repEig.dcm']);
            dicomwrite(uint8(255*LA_cfg_hsv),[foutputdir,name,'_3-9_drLA-cfg2-',num2str(nr),'repEig_hsvColoring.dcm']);
            dicomwrite(PRRc,[foutputdir,name,'_3-10_PhR-cfg2-',num2str(nr),'repEig.dcm']);
        end
        end

    end% save results
    fprintf('文件 %s 处理完成!\n', display_name);
    fprintf('输出目录: %s\n', foutputdir);
    fprintf('\n');
end

end

%%
return;


%%
% cal QUV
function [S0,S1,S2,S3] = cumulativeQUV(IMG_ch1,IMG_ch2)
axis=angle(IMG_ch2.*conj(IMG_ch1));
S0=abs(IMG_ch1).^2+abs(IMG_ch2).^2;
S1=abs(IMG_ch1).^2-abs(IMG_ch2).^2;
S2=2.*abs(IMG_ch1).*abs(IMG_ch2).*cos(axis);
S3=2.*abs(IMG_ch1).*abs(IMG_ch2).*sin(-axis);
end

function [OAC] = calOAC(linFrame)
OAC=linFrame*0;
Nois_M=mean(mean(linFrame(size(linFrame,1)-2:size(linFrame,1),:)));
if Nois_M < 1, return; end
Nois_D=std(std(linFrame(size(linFrame,1)-2:size(linFrame,1),:)));
linFrame=linFrame-Nois_M+Nois_D; %去除噪声地板
down=linFrame<0;linFrame(down)=1; %去除负数部分
tail_signal=mean(linFrame(size(linFrame,1)-5:size(linFrame,1),:)); %估计剩余光强
tail_signal = medfilt1(tail_signal,15);
tail_signal = smooth(tail_signal, 15)';

for z=1:size(linFrame,1)
%     OAC(z,:)=linFrame(z,:)./(2*0.005*sum(linFrame(z+1:size(linFrame,1),:))+tail_signal);
    OAC(z,:)=linFrame(z,:)./(2*0.0086*sum(linFrame(z+1:size(linFrame,1),:))+tail_signal);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% surface segmentation %%%%%%%%%%%%%%%%%%%%%%%%%
% 使用OAC数据检测样品，首先计算梯度，然后检测第一个峰值
% OAC: 输入的OAC数据
% surf_threshold: 检测大于这个阈值的第一个峰值

function [test_seg_top]=surf_seg(OAC,surf_threshold)

test_seg_top = ones(1,size(OAC,2));
if sum(OAC(:)) < 1,return;end
h1=[-0.2;-0.2;-0.2;-0.2;-0.2;1]; %上表面梯度模版
OAC_G=(filter2(h1,OAC)>surf_threshold); %上表面模板滤波
se = strel('disk', 2); % 选择一个半径为5的圆形结构元素
OAC_G = imdilate(OAC_G, se);
OAC_G = imerode(OAC_G, se);
OAC_G = double(bwareaopen(OAC_G, 20));
%figure,imshow(OAC_G,[0 1]);
Temp_locs=1;
for j=1:size(OAC,2)
        if sum(OAC_G(1:3,j),'all') > 1, ILM(j) = 1; Temp_locs=1;continue; end
[~,locs] = findpeaks(OAC_G(:,j),'minpeakheight',0.2);% mph_ILM 设定峰值的最小高度
        if isempty(locs) 
            locs(1)=Temp_locs;
        end
        ILM(j)=locs(1)+1; 
        Temp_locs=locs(1);
end
ILM = medfilt1(ILM,11);
test_seg_top = round(smooth(ILM', 15))+1;

end


function [LA_c,PhR_c,cumLA] = calLAPhR(IMG_ch1,IMG_ch2,test_seg_top,dopu_splitSpec_M,kRL,kRU,h1,h2,Avnum,wovWinF)
[nZ,nX] = size(IMG_ch1);
LA_c = zeros(nZ-Avnum,nX,3);PhR_c = zeros(nZ-Avnum,nX);
cumLA = zeros(nZ-Avnum,nX,3);
[ES0,ES1,ES2,ES3] = cumulativeQUV(IMG_ch1,IMG_ch2);
if sum(ES0(:)) < 5, return;end
EQm=ES1./ES0;
EUm=ES2./ES0;
EVm=ES3./ES0;

Stru_E = zeros(nZ,nX);

if wovWinF
    EQmm=imfilter(EQm,h1,'replicate');
    EUmm=imfilter(EUm,h1,'replicate');
    EVmm=imfilter(EVm,h1,'replicate');
else
    [EQmm] = vWinAvgFiltOpt_2_1(EQm,dopu_splitSpec_M,kRL,kRU);
    [EUmm] = vWinAvgFiltOpt_2_1(EUm,dopu_splitSpec_M,kRL,kRU);
    [EVmm] = vWinAvgFiltOpt_2_1(EVm,dopu_splitSpec_M,kRL,kRU);
end
[LA_c,PhR_c,cumLA] = FreeSpace_PSOCT_3_DDG_rmBG_7(EQmm,EUmm,EVmm,Stru_E,test_seg_top,h1,h2,Avnum);
end

%% drLA with data do zero-avg loacally
function [LA_c,PhR_c,cumLA] = calLAPhRcfg2(IMG_ch1,IMG_ch2,test_seg_top,dopu_splitSpec_M,kRL,kRU,h1,h2,Avnum,wovWinF)
[nZ,nX] = size(IMG_ch1);
LA_c = zeros(nZ-Avnum,nX,3);PhR_c = zeros(nZ-Avnum,nX);
cumLA = zeros(nZ-Avnum,nX,3);
[ES0,ES1,ES2,ES3] = cumulativeQUV(IMG_ch1,IMG_ch2);
if sum(ES0(:)) < 5, return;end
EQm=ES1./ES0;
EUm=ES2./ES0;
EVm=ES3./ES0;

Stru_E = zeros(nZ,nX);

if wovWinF
    EQmm=imfilter(EQm,h1,'replicate');
    EUmm=imfilter(EUm,h1,'replicate');
    EVmm=imfilter(EVm,h1,'replicate');
else
    [EQmm] = vWinAvgFiltOpt_2_1(EQm,dopu_splitSpec_M,kRL,kRU);
    [EUmm] = vWinAvgFiltOpt_2_1(EUm,dopu_splitSpec_M,kRL,kRU);
    [EVmm] = vWinAvgFiltOpt_2_1(EVm,dopu_splitSpec_M,kRL,kRU);
end
% DDG to cal drLA without BG modulation
[LA_c,PhR_c,cumLA] = FreeSpace_PSOCT_3_DDG_rmBG_7(EQmm,EUmm,EVmm,Stru_E,test_seg_top,h1,h2,Avnum);

end
%
function [LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = calLAPhRALL(IMG_ch1,IMG_ch2,test_seg_top,dopu_splitSpec_M,kRL,kRU,h1,h2,Avnum,wovWinF)

if nargin < 10, wovWinF = 0;end
[nZ,nX] = size(IMG_ch1);
LA = zeros(nZ-Avnum,nX,3);PhR = zeros(nZ-Avnum,nX);cumLA = LA;
LA_raw = zeros(nZ-Avnum,nX,3);PhR_raw = zeros(nZ-Avnum,nX);cumLA_raw = LA;
[ES0,ES1,ES2,ES3] = cumulativeQUV(IMG_ch1,IMG_ch2);
if sum(ES0(:)) < 5, return;end
EQm=ES1./ES0;
EUm=ES2./ES0;
EVm=ES3./ES0;
% QUV_E(:,:,:)=cat(3,EQm,EUm,EVm);
% Stru_E=20*log10(S0);
Stru_E = zeros(nZ,nX);
if wovWinF
EQmm=imfilter(EQm,h1,'replicate'); 
EUmm=imfilter(EUm,h1,'replicate'); 
EVmm=imfilter(EVm,h1,'replicate');   
else
[EQmm] = vWinAvgFiltOpt_2_1(EQm,dopu_splitSpec_M,kRL,kRU);
[EUmm] = vWinAvgFiltOpt_2_1(EUm,dopu_splitSpec_M,kRL,kRU);
[EVmm] = vWinAvgFiltOpt_2_1(EVm,dopu_splitSpec_M,kRL,kRU);
end

[LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = FreeSpace_PSOCT_3_DDG_rmBG_7(EQmm,EUmm,EVmm,Stru_E,test_seg_top,h1,h2,Avnum);
end
%% flat dcm results
function [LAs_flat] = volFlatten(LAs,topLines)
[nZ,nX,nCh,nY] = size(LAs);
LAs_flat = zeros(nZ,nX,nCh,nY);
for i = 1:nX
   for j= 1:nY
       LAs_flat(1:nZ - topLines(i,j)+1,i,:,j) = ...
           LAs(topLines(i,j):nZ,i,:,j);
   end
end
end

%% QU recoloring
function [hsvLA] = quColoring(fLAs,cmapshift)
    showimg = 0;
    cusmapRg = hsv(512);
%     cusmapRg = cusColormap(0);
    if nargin > 1, cusmapRg = circshift(cusmapRg,cmapshift); end
%     cusmapRg = [hsv(256);hsv(256)];%
    thetaRg = linspace(-pi,pi,256*2);
    thetas = atan2(fLAs(:,:,2),fLAs(:,:,1));

    
    thpf = polyfit(thetaRg,1:256*2,1);
    thetasInds = round(polyval(thpf,thetas));
    colorInds = cusmapRg(thetasInds,:);
    hsvLA = reshape(colorInds,[size(thetasInds),3]);
    if showimg
        figure;imshow(hsvLA,[]);
    end
end