%%% rPSOCT_process_single_file.m
%%% 1. frame based code to cal 3D LAs and Phase[two cfg: cfg1== LA, cfg2==Phase]
%%% 1. use doup to apply a variable gaussian fitler to Q U V before cal LA
%%%   and phase
%%% 2. speed the cal process, params.range.setZrg = 100; params.range.setZrg = 0 for process all
%%% ______________________________20240913_______________________________________________
%%% 3. optional to toggle vWinF
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

% 设置输入输出路径
data_path   = 'D:\1-Liu Jian\yongxin.wang\tmp\';
output_base = 'D:\1-Liu Jian\yongxin.wang\tmp1';
if ~exist(data_path, 'dir')
    error(['数据路径不存在: ' data_path]);
end

% 创建输出路径
if ~exist(output_base, 'dir')
    fprintf(['输出路径不存在: ' output_base]);
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
            % 记录开始时间
            tic;
            
            % 调用单文件处理函数
            rPSOCT_process_single_file(full_path, output_base);
            
            % 计算并显示处理时间
            proc_time = toc;
            fprintf('文件 %s 处理成功, 耗时: %.2f 秒 (%.2f 分钟)\n', oct_files(i).name, proc_time, proc_time/60);
        catch ME
            fprintf('处理文件 %s 时出错: %s\n', oct_files(i).name, ME.message);
            % 继续处理下一个文件
        end
    end
    
    fprintf('\n所有文件处理完成!\n');
    
end


function rPSOCT_process_single_file(varargin)
    clc;
    %% (1) 获取输入参数
    if nargin < 1
        error('请提供输入文件路径作为参数');
    end
    
    input_file_path = varargin{1};
    
    % 如果提供了输出路径，则使用提供的路径，否则使用默认路径
    if nargin >= 2
        output_base = varargin{2};
    else
        output_base = 'D:\1-Liu Jian\yongxin.wang\PSOCT\2025-9-19\ssdopu-kRL_16-kRU_23';
    end
    
    % 检查输入文件是否存在
    if ~exist(input_file_path, 'file')
        error(['输入文件不存在: ' input_file_path]);
    end
    
    % 创建输出目录
    if ~exist(output_base,'dir'), mkdir(output_base), end
    
    % 记录函数开始时间
    file_start_time = tic;
    
    %% (2) 初始化参数
    
    % 启动并行池
    if isempty(gcp('nocreate'))
        parpool('local', 48);
    end
    
    params = config_params();
    
    
    % 直接用params.xxx结构体成员，无需单独赋值
    distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing
    
    %% （4）处理单个文件
    filename = input_file_path;
    [filepath, name, ext] = fileparts(filename);
    display_name = [name ext];
    
    fprintf('\n==================================================\n');
    fprintf('正在处理文件: %s\n', display_name);
    fprintf('文件路径: %s\n', filename);
    fprintf('==================================================\n');
    
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
    windex=1 : winL / 2 : Blength;
    
    if params.processing.useref==1 %use the first 50k a-lines to calc ref
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
    elseif params.processing.useref==-1
        Ref_ch1=0;
        Ref_ch2=0;
    else
        Ref_ch1=repmat(Bck1,[1,nX]);
        Ref_ch2=repmat(Bck2,[1,nX]);
    end
    
    
    %(b) Calculates the Cumulative Stokes parameters I,Q,U,V
    for nr = nR 
        if params.range.setZrg
            foutputdir = fullfile(output_base,[name,'.rep',num2str(nr),'fastZ_testDDG',num2str(params.range.setZrg),'\']);
        else
            foutputdir = fullfile(output_base,[name,'.rep',num2str(nr),'DDG_corr3_wobg_test\']);
        end
        
        if ~exist(foutputdir,'dir')
            mkdir(foutputdir); 
        end
        
        czrg = 1:320;% set z range
        topLines = ones(nX,nY);
        
        if params.processing.hasSeg
            % 尝试从与输入文件相同的文件夹加载 .mat（避免依赖当前工作目录）
            matpath = fullfile(filepath, [name, '.mat']);
            if exist(matpath, 'file')
            load(matpath);
            fprintf('成功加载边界mat文件: %s\n', matpath);
            topLines = round(topLines); 
            topLines(topLines<=1) = 1;
        else
            warning('在预计位置未找到 %s，将使用自动分割代替。', matpath);
        end
    end
    
    if params.range.setZrg
        czrg = czrg(1:round(params.range.setZrg));
    end
    nZcrop = numel(czrg);nZ = nZcrop;
    % (c) creat array to store results
    LA_c_cfg1_avg = zeros(nZcrop-params.polarization.Avnum,nX,3,nY);
    PhR_c_cfg1_avg = zeros(nZcrop-params.polarization.Avnum,nX,nY);
    cumLA_cfg1_avg = zeros(nZcrop-params.polarization.Avnum,nX,3,nY);
    
    LA_Ms_cfg1_rmBG = zeros(nZcrop-params.polarization.Avnum,nX,3,nY);
    PhR_Ms_cfg1_rmBG = zeros(nZcrop-params.polarization.Avnum,nX,3,nY);
    cumLA_Ms_cfg1_rmBG = zeros(nZcrop-params.polarization.Avnum,nX,3,nY);
    %% cfg2
    Smap_avg = zeros(numel(czrg),nX,3,nY);Smap_rep1 = zeros(numel(czrg),nX,3,nY);
    Strus = zeros(numel(czrg),nX,nY);Stru_OAC = Strus;
    dopu_splitSpectrum = zeros(numel(czrg),nX,nY);
    fprintf('开始处理 %d 个B-Scan...\n', nY);
    parfor iY=1:nY % Bscan number (parfor不支持指定步长)
        
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
        if params.processing.do_PhComp==1
            Bd1 =real(hilbert(Bs1).*phV);  % dispersion correction
            Bd2 =real(hilbert(Bs2).*phV);
        else
            Bd1=Bs1;
            Bd2=Bs2;
        end
        %% split spectrum
        if ~params.dopu.do_ssdopu % if not do adaptive gaussian filter dopu_ss =1
            dopu_ss = 1;
        else
            % creat array to store split-spectrum complex (FFT): Z*X*nR*nWin
            Bimg1 = zeros(SPL,nX,nr,nWin);
            Bimg2 = Bimg1;
            S0 = zeros(nZcrop,nX,nr,nWin); 
            S1 = S0;
            S2 = S0;
            S3 = S0;
            for iR = 1:nr
                for iL=1:nWin
                    % extract data fragments from two different channel and
                    % apply a giassian filter before performing FFT
                    iBd1=Bd1(windex(iL):windex(iL)+winL-1,:,iR).*winG;
                    iBd2=Bd2(windex(iL):windex(iL)+winL-1,:,iR).*winG;
                    Bimg1(:,:,iR,iL)=fft(iBd1,SPL,1);
                    Bimg2(:,:,iR,iL)=fft(iBd2,SPL,1);
                end
            end
            IMGs_ch1=Bimg1(czrg,:,:,:);IMGs_ch2=Bimg2(czrg,:,:,:);
            % calculate the QUV from the two-channel complex
            for iR = 1:nr
                for iL = 1: nWin
                    [S0(:,:,iR,iL),S1(:,:,iR,iL),S2(:,:,iR,iL),S3(:,:,iR,iL)] = ...
                        cumulativeQUV(IMGs_ch1(:,:,iR,iL),IMGs_ch2(:,:,iR,iL));
                end
            end
            % average QUV accross split-spectrum and calculate dopu
            dopu_ss = sqrt(mean(S1./S0,4).^2+mean(S2./S0,4).^2+mean(S3./S0,4).^2);%% ===>nZcrop,nX,nr
            dopu_splitSpectrum(:,:,iY) = mean(dopu_ss,3);
        end
        %% whole spectrum nZ*nX*nR==> fft(complex)
        Bimg1_wholeStr=fft(Bd1.*winG_whole,SPL,1);Bimg2_wholeStr=fft(Bd2.*winG_whole,SPL,1);
        IMG1_wholeStr = Bimg1_wholeStr(czrg,:,:);IMG2_wholeStr = Bimg2_wholeStr(czrg,:,:);
        
        %% Struc Stokes, and OAC
        [wS0,wS1,wS2,wS3] = cumulativeQUV(IMG1_wholeStr,IMG2_wholeStr);
        wQ = wS1./wS0;wU = wS2./wS0; wV = wS3./wS0;
        strLin = mean(wS0,3);
        Strus(:,:,iY) = 20*log10(strLin);
        Smap_avg(:,:,:,iY) = cat(3,mean(wQ,3),mean(wU,3),mean(wV,3));
        Smap_rep1(:,:,:,iY) = cat(3,wQ(:,:,1),wU(:,:,1),wV(:,:,1));
        if ~params.processing.hasSeg, strOAC = calOAC(strLin); topLines(:,iY)=surf_seg(strOAC,0.25)+2; end
        
        %% drLA and drPhR
        dopu_splitSpec_M = squeeze(mean(dopu_ss,3)); %% dopu across the different repeat(nR)
        if params.mode.do_cfg1
            if params.dopu.do_avg %% avg -- cfg1
                IMG_ch1 = squeeze(mean(IMG1_wholeStr,3));IMG_ch2 = squeeze(mean(IMG2_wholeStr,3));
                [LA_c_cfg1_avg(:,:,:,iY),PhR_c_cfg1_avg(:,:,iY),cumLA_cfg1_avg(:,:,:,iY),...
                    LA_Ms_cfg1_rmBG(:,:,:,iY),PhR_Ms_cfg1_rmBG(:,:,iY),cumLA_Ms_cfg1_rmBG(:,:,:,iY)] = ...
                    calLAPhRALL(IMG_ch1,IMG_ch2,topLines(:,iY),dopu_splitSpec_M,params.polarization.kRL_cfg1,params.polarization.kRU_cfg1,params.filters.h1,params.filters.h2,params.polarization.Avnum,params.mode.wovWinF);
            end
        end
        fclose(fid);
    end % parfor iY
    fclose all;
    %% save results: strus(flow),stokes,oac
    if params.tiff.saveDicom
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
        
        % 对DOPU应用阈值处理：使用边界信息，边界以上的区域视为噪声并设为纯黑；
        % 边界到向下的前100个坐标（包含表面位置）不被过滤，100像素以下开始应用阈值过滤
        dopu_threshold = 0.5;

        % 创建DOPU过滤结果：初始化为原始值
        dopu_filtered = dopu_splitSpectrum;

        [nZ_dopu, nX_dopu, nY_dopu] = size(dopu_splitSpectrum);

        for iY = 1:nY_dopu
            for iX = 1:nX_dopu
                % 获取该A-line的表面位置
                surface_pos = topLines(iX, iY);

                % 边界以上的区域（噪声区域）设为纯黑
                if surface_pos > 1
                    noise_range = 1:surface_pos-1;
                    dopu_filtered(noise_range, iX, iY) = 0;  % 纯黑
                end

                % 边界以下的前100个像素（包含表面）保持原值，不应用过滤
                if surface_pos <= nZ_dopu
                    keep_start = surface_pos;
                    keep_end = min(nZ_dopu, surface_pos + 100 - 1); % 前100个坐标
                    % 从 keep_end+1 开始才应用阈值过滤
                    start_filter = keep_end + 1;
                    if start_filter <= nZ_dopu
                        filter_range = start_filter:nZ_dopu;
                        dopu_filtered(filter_range, iX, iY) = dopu_splitSpectrum(filter_range, iX, iY) .* (dopu_splitSpectrum(filter_range, iX, iY) >= dopu_threshold);
                    end
                end
            end
        end

        dicomwrite(uint8(255*(permute(dopu_splitSpectrum,[1 2 4 3]))),[foutputdir,name,'_1-4_dopu_SS.dcm']);
        dicomwrite(uint8(255*(permute(dopu_filtered,[1 2 4 3]))),[foutputdir,name,'_1-5_dopu_SS_filtered.dcm']);
        
        
        if ~params.processing.hasSeg
            save([foutputdir,name, '.mat'],'topLines','czrg');
        end
        save([foutputdir,'topLines.mat'],'topLines','czrg');
        rotAngle = 440;
        if params.mode.do_cfg1
            PRRrg = [0 0.5];
            writematrix(PRRrg,[foutputdir,name,'_2-0_PhRRg.txt']);
            if params.dopu.do_avg
                % 创建DOPU掩码：DOPU > 0.6的位置为1，否则为0
                dopu_mask = dopu_splitSpectrum > dopu_threshold;

                % 裁剪DOPU掩码以匹配偏振参数的维度（减去平均层数）
                dopu_mask_cropped = dopu_mask(1:end-params.polarization.Avnum, :, :);

                % 修改掩码：使用边界信息，边界以上的区域视为噪声并设为纯黑，边界以下应用阈值过滤
                [nZ_cropped, nX, nY] = size(dopu_mask_cropped);

                % 创建新的掩码：边界以上的噪声区域设为0（纯黑），边界以下应用DOPU阈值
                dopu_mask_bottom_only = zeros(size(dopu_mask_cropped));

                for iY = 1:nY
                    for iX = 1:nX
                        % 获取该A-line的表面位置
                        surface_pos = topLines(iX, iY);

                        % 计算相对于裁剪后数据的表面位置
                        surface_pos_cropped = max(1, surface_pos - params.polarization.Avnum);

                        % 边界以上的区域（噪声区域）设为纯黑
                        if surface_pos_cropped > 1
                            noise_range = 1:surface_pos_cropped-1;
                            dopu_mask_bottom_only(noise_range, iX, iY) = 0;  % 纯黑
                        end

                                % 边界以下的前100个像素（包含表面）保持不被过滤（掩码为1），之后才使用原始掩码
                                if surface_pos_cropped <= nZ_cropped
                                    keep_start_c = surface_pos_cropped;
                                    keep_end_c = min(nZ_cropped, surface_pos_cropped + 100 - 1); % 包含表面，共100个坐标
                                    % 将前100像素设为1（不滤）
                                    dopu_mask_bottom_only(keep_start_c:keep_end_c, iX, iY) = 1;
                                    % 对更深的像素应用原始掩码
                                    start_apply_c = keep_end_c + 1;
                                    if start_apply_c <= nZ_cropped
                                        deeper_range = start_apply_c:nZ_cropped;
                                        dopu_mask_bottom_only(deeper_range, iX, iY) = dopu_mask_cropped(deeper_range, iX, iY);
                                    end
                                end
                    end
                end
                
                for iY =1:nY
                    % 原始PhR和LA处理
                    PRRc(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_avg(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_cfg_hsv(:,:,:,iY)] = quColoring(cumLA_cfg1_avg(:,:,:,iY),rotAngle);
                    [LA_cfg_hsv(:,:,:,iY)] = quColoring(LA_c_cfg1_avg(:,:,:,iY),rotAngle);
                    
                    PRRc_rmBG(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_Ms_cfg1_rmBG(:,:,iY),PRRrg)*256),parula(256))*256);
                    [cumLA_Ms_cfg1_rmBG_hsv(:,:,:,iY)] = quColoring(cumLA_Ms_cfg1_rmBG(:,:,:,iY),rotAngle);
                    [LA_Ms_cfg1_rmBG_hsv(:,:,:,iY)] = quColoring(LA_Ms_cfg1_rmBG(:,:,:,iY),rotAngle);
                    
                    % 应用DOPU阈值过滤：只显示DOPU > 0.6的区域
                    % 过滤后的PhR和LA（低于阈值的区域设为0）
                    PhR_c_cfg1_avg_filtered = PhR_c_cfg1_avg(:,:,iY) .* dopu_mask_bottom_only(:,:,iY);
                    PhR_Ms_cfg1_rmBG_filtered = PhR_Ms_cfg1_rmBG(:,:,iY) .* dopu_mask_bottom_only(:,:,iY);
                    
                    cumLA_cfg1_avg_filtered = cumLA_cfg1_avg(:,:,:,iY);
                    cumLA_cfg1_avg_filtered(repmat(~dopu_mask_bottom_only(:,:,iY), [1 1 3])) = 0;
                    
                    LA_c_cfg1_avg_filtered = LA_c_cfg1_avg(:,:,:,iY);
                    LA_c_cfg1_avg_filtered(repmat(~dopu_mask_bottom_only(:,:,iY), [1 1 3])) = 0;
                    
                    cumLA_Ms_cfg1_rmBG_filtered = cumLA_Ms_cfg1_rmBG(:,:,:,iY);
                    cumLA_Ms_cfg1_rmBG_filtered(repmat(~dopu_mask_bottom_only(:,:,iY), [1 1 3])) = 0;
                    
                    LA_Ms_cfg1_rmBG_filtered = LA_Ms_cfg1_rmBG(:,:,:,iY);
                    LA_Ms_cfg1_rmBG_filtered(repmat(~dopu_mask_bottom_only(:,:,iY), [1 1 3])) = 0;
                    
                    % 生成过滤后的彩色图像
                    PRRc_filtered(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_c_cfg1_avg_filtered,PRRrg)*256),parula(256))*256);
                    [cumLA_cfg_hsv_filtered(:,:,:,iY)] = quColoring(cumLA_cfg1_avg_filtered,rotAngle);
                    [LA_cfg_hsv_filtered(:,:,:,iY)] = quColoring(LA_c_cfg1_avg_filtered,rotAngle);
                    
                    PRRc_rmBG_filtered(:,:,:,iY) = uint8(ind2rgb(uint8(mat2gray(PhR_Ms_cfg1_rmBG_filtered,PRRrg)*256),parula(256))*256);
                    [cumLA_Ms_cfg1_rmBG_hsv_filtered(:,:,:,iY)] = quColoring(cumLA_Ms_cfg1_rmBG_filtered,rotAngle);
                    [LA_Ms_cfg1_rmBG_hsv_filtered(:,:,:,iY)] = quColoring(LA_Ms_cfg1_rmBG_filtered,rotAngle);
                    
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
                
                % 输出过滤后的DCM文件（2-11到2-20系列）
                dicomwrite(uint8(255*(cumLA_cfg1_avg_filtered/2+0.5)),[foutputdir,name,'_2-11_cumLA_filtered-cfg1-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*cumLA_cfg_hsv_filtered),[foutputdir,name,'_2-12_cumLA_filtered-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_c_cfg1_avg_filtered/2+0.5)),[foutputdir,name,'_2-13_drLA_filtered-cfg1-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*LA_cfg_hsv_filtered),[foutputdir,name,'_2-14_drLA_filtered-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(PRRc_filtered,[foutputdir,name,'_2-15_PhR_filtered-cfg1-',num2str(nr),'repAvg.dcm']);
                
                dicomwrite(uint8(255*(cumLA_Ms_cfg1_rmBG_filtered/2+0.5)),[foutputdir,name,'_2-16_cumLA_rmBG_filtered-cfg1-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*cumLA_Ms_cfg1_rmBG_hsv_filtered),[foutputdir,name,'_2-17_cumLA_rmBG_filtered-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(uint8(255*(LA_Ms_cfg1_rmBG_filtered/2+0.5)),[foutputdir,name,'_2-18_drLA_rmBG_filtered-cfg1-',num2str(nr),'repAvg.dcm']);
                dicomwrite(uint8(255*LA_Ms_cfg1_rmBG_hsv_filtered),[foutputdir,name,'_2-19_drLA_rmBG_filtered-cfg1-',num2str(nr),'repAvg_hsvColoring.dcm']);
                dicomwrite(PRRc_rmBG_filtered,[foutputdir,name,'_2-20_PhR_rmBG_filtered-cfg1-',num2str(nr),'repAvg.dcm']);
            end
        end
    end
    %% DCM to TIFF
    if params.tiff.saveDicom && params.tiff.make_tiff
        fprintf('\n开始从DCM提取TIFF\n');
        try
            % 调用本地DCM到TIFF转换函数
            convert_dcm_to_tiff_local(foutputdir, params.tiff.tiff_frame, name);
        catch ME
            fprintf('DCM到TIFF转换出错: %s\n', ME.message);
            % 转换出错不影响主流程，继续执行
        end
    end
    
    fprintf('文件 %s 处理完成!\n', display_name);
    fprintf('输出目录: %s\n', foutputdir);
    
    % 显示处理时间统计
    file_proc_time = toc(file_start_time);
    fprintf('处理时间: %.2f 秒 (%.2f 分钟)\n', file_proc_time, file_proc_time/60);
    fprintf('\n');
    
    % 清理并行池（在函数结束时）
    delete(gcp('nocreate'));
    end
return;


%% ====================================================================================
% 函数名: cumulativeQUV
% 功能: 从两通道复数OCT图像计算Stokes参数(S0, S1, S2, S3)
% 输入参数:
%   IMG_ch1 - 通道1的复数OCT信号矩阵 [Z×X或Z×X×Rep]
%   IMG_ch2 - 通道2的复数OCT信号矩阵 [Z×X或Z×X×Rep]
% 输出参数:
%   S0 - Stokes参数S0 (总光强度): |E1|² + |E2|²
%   S1 - Stokes参数S1 (水平/垂直偏振差): |E1|² - |E2|²
%   S2 - Stokes参数S2 (±45°偏振差): 2|E1||E2|cos(θ)
%   S3 - Stokes参数S3 (左/右圆偏振差): 2|E1||E2|sin(-θ)
% 说明:
%   - Stokes参数是描述偏振光的完整参数集
%   - θ为两通道信号的相位差
%   - 用于后续DOPU(偏振均匀度)和偏振分析计算
% ====================================================================================
function [S0,S1,S2,S3] = cumulativeQUV(IMG_ch1,IMG_ch2)
    % 计算两通道信号的相位差
    axis = angle(IMG_ch2.*conj(IMG_ch1));

    % 计算四个Stokes参数
    S0 = abs(IMG_ch1).^2 + abs(IMG_ch2).^2;                     % 总光强度
    S1 = abs(IMG_ch1).^2 - abs(IMG_ch2).^2;                     % 水平-垂直偏振分量差
    S2 = 2.*abs(IMG_ch1).*abs(IMG_ch2).*cos(axis);              % +45°与-45°偏振分量差
    S3 = 2.*abs(IMG_ch1).*abs(IMG_ch2).*sin(-axis);             % 右旋与左旋圆偏振分量差
end

%% ====================================================================================
% 函数名: calOAC
% 功能: 计算光学衰减系数(Optical Attenuation Coefficient)
% 输入参数:
%   linFrame - 线性OCT信号强度矩阵 [Z×X]，通常为|FFT(OCT信号)|²
% 输出参数:
%   OAC - 光学衰减系数矩阵 [Z×X]，用于组织边界检测和结构分析
% 算法原理:
%   - 基于Beer-Lambert定律: I(z) = I₀ exp(-2μz)
%   - μ为衰减系数，反映组织的散射和吸收特性
%   - 通过当前信号强度与剩余信号总和的比值估算局部衰减
% 应用: 主要用于组织表面检测、分层分析和结构图像增强
% ====================================================================================
function [OAC] = calOAC(linFrame)
    % 初始化输出矩阵
    OAC = linFrame * 0;

    % 估算噪声底限(使用最后3层的均值)
    Nois_M = mean(mean(linFrame(size(linFrame,1)-2:size(linFrame,1),:)));
    if Nois_M < 1, return; end  % 信号过弱则直接返回

    % 计算噪声标准差
    Nois_D = std(std(linFrame(size(linFrame,1)-2:size(linFrame,1),:)));

    % 去除噪声底限，添加噪声标准差作为最小信号阈值
    linFrame = linFrame - Nois_M + Nois_D;

    % 将负值设为1，避免对数运算错误
    down = linFrame < 0;
    linFrame(down) = 1;

    % 估算尾部残余光强(最后6层的平均值)
    tail_signal = mean(linFrame(size(linFrame,1)-5:size(linFrame,1),:));
    tail_signal = medfilt1(tail_signal, 15);    % 中值滤波去除突变
    tail_signal = smooth(tail_signal, 15)';      % 平滑滤波

    % 计算每一深度层的光学衰减系数
    % OAC = I(z) / (2α∑I(z':z') + 残余光强)
    % 其中α=0.0086为经验衰减参数
    for z = 1:size(linFrame,1)
        OAC(z,:) = linFrame(z,:) ./ (2*0.0086*sum(linFrame(z+1:size(linFrame,1),:)) + tail_signal);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% 组织表面分割算法 %%%%%%%%%%%%%%%%%%%%%%%%%
%% ====================================================================================
% 函数名: surf_seg
% 功能: 基于OAC数据进行组织表面检测和分割
% 输入参数:
%   OAC - 光学衰减系数矩阵 [Z×X]
%   surf_threshold - 表面检测阈值(通常为0.25)，用于区分组织与背景
% 输出参数:
%   test_seg_top - 每个A-line的组织表面位置向量 [1×X]
% 算法流程:
%   1. 使用梯度模板检测表面边界
%   2. 形态学处理去除噪声
%   3. 峰值检测确定表面位置
%   4. 中值滤波和平滑处理获得连续表面
% 应用: 用于C-scan展平、分层分析和深度校正
% ====================================================================================

function [test_seg_top] = surf_seg(OAC, surf_threshold)
    % 初始化表面位置向量(默认值为1)
    test_seg_top = ones(1, size(OAC, 2));

    % 检查输入数据有效性
    if sum(OAC(:)) < 1, return; end

    % 定义上表面检测的梯度模板(检测从低到高的跳变)
    params.filters.h1 = [-0.2; -0.2; -0.2; -0.2; -0.2; 1]; % 5个负权重 + 1个正权重

    % 使用梯度模板进行卷积，检测超过阈值的边界
    OAC_G = (filter2(params.filters.h1, OAC) > surf_threshold);

    % 形态学处理：先膨胀后腐蚀(闭运算)，连接断裂边界
    se = strel('disk', 2);          % 创建半径为2的圆形结构元素
    OAC_G = imdilate(OAC_G, se);    % 膨胀操作
    OAC_G = imerode(OAC_G, se);     % 腐蚀操作

    % 去除小于20像素的连通区域(噪声去除)
    OAC_G = double(bwareaopen(OAC_G, 20));

    % 备用显示: figure, imshow(OAC_G, [0 1]);

    % 逐列检测表面位置
    Temp_locs = 1;  % 临时位置变量，用于处理检测失败的情况

    for j = 1:size(OAC, 2)
        % 检查前3行是否有信号(防止表面过于靠近顶部)
        if sum(OAC_G(1:3, j), 'all') > 1
            ILM(j) = 1;
            Temp_locs = 1;
            continue;
        end
        
        % 使用峰值检测找到第一个显著峰值
        [~, locs] = findpeaks(OAC_G(:, j), 'minpeakheight', 0.2);
        
        % 如果没有检测到峰值，使用前一列的位置
        if isempty(locs)
            locs(1) = Temp_locs;
        end
        
        % 记录检测到的表面位置(+1偏移避免边界问题)
        ILM(j) = locs(1) + 1;
        Temp_locs = locs(1);
    end

    % 中值滤波去除突变点(窗口大小11)
    ILM = medfilt1(ILM, 11);

    % 平滑滤波获得连续表面(窗口大小15)，并添加安全偏移+1
    test_seg_top = round(smooth(ILM', 15)) + 1;
end

%% ====================================================================================
% 函数名: calLAPhRALL (完整版偏振参数计算)
% 功能: 计算完整的偏振参数，包括处理前和处理后的结果
% 输入参数: (同前述函数)
%   IMG_ch1, IMG_ch2 - 双通道复数OCT信号
%   test_seg_top - 组织表面位置
%   dopu_splitSpec_M - 分光谱DOPU矩阵
%   kRL, kRU - 滤波核范围参数
%   params.filters.h1, params.filters.h2 - 多尺度高斯滤波核
%   params.polarization.Avnum - 平均层数
%   params.mode.wovWinF - 滤波模式(可选参数，默认为0)
% 输出参数:
%   LA - 处理后的局部双折射 [(Z-params.polarization.Avnum)×X×3]
%   PhR - 处理后的相位延迟 [(Z-params.polarization.Avnum)×X]
%   cumLA - 处理后的累积双折射 [(Z-params.polarization.Avnum)×X×3]
%   LA_raw - 原始局部双折射(去背景前)
%   PhR_raw - 原始相位延迟(去背景前)
%   cumLA_raw - 原始累积双折射(去背景前)
% 特点:
%   - 提供完整的处理流程，包括原始和处理后的结果
%   - 支持两种滤波模式的完整实现
%   - 用于深入的偏振分析和算法验证
% ====================================================================================
function [LA,PhR,cumLA,LA_raw,PhR_raw,cumLA_raw] = calLAPhRALL(IMG_ch1,IMG_ch2,test_seg_top,dopu_splitSpec_M,kRL,kRU,h1,h2,Avnum,wovWinF)
    % 处理可选参数
    if nargin < 10, wovWinF = 0; end  % 直接设置wovWinF而不是params.mode.wovWinF

    % 获取数据维度
    [nZ, nX] = size(IMG_ch1);

    % 初始化处理后的输出矩阵
    LA = zeros(nZ-Avnum, nX, 3);       % 使用传入的Avnum参数
    PhR = zeros(nZ-Avnum, nX);         % 使用传入的Avnum参数
    cumLA = LA;                        % 处理后累积双折射

    % 初始化原始输出矩阵
    LA_raw = zeros(nZ-Avnum, nX, 3);   % 使用传入的Avnum参数
    PhR_raw = zeros(nZ-Avnum, nX);     % 使用传入的Avnum参数
    cumLA_raw = LA;                    % 原始累积双折射

    % 计算Stokes参数
    [ES0, ES1, ES2, ES3] = cumulativeQUV(IMG_ch1, IMG_ch2);

    % 信号强度检查
    if sum(ES0(:)) < 5, return; end

    % 计算归一化Stokes分量
    EQm = ES1 ./ ES0;   % Q分量(水平-垂直偏振差)
    EUm = ES2 ./ ES0;   % U分量(±45°偏振差)
    EVm = ES3 ./ ES0;   % V分量(左右旋圆偏振差)

    % 备用代码: QUV_E(:,:,:) = cat(3, EQm, EUm, EVm);
    % 备用代码: Stru_E = 20*log10(S0);

    % 初始化结构矩阵
    Stru_E = zeros(nZ, nX);

    % 根据滤波模式进行预处理
    if wovWinF == 1
        % 模式1: 固定高斯滤波
        EQmm = imfilter(EQm, h1, 'replicate');
        EUmm = imfilter(EUm, h1, 'replicate');
        EVmm = imfilter(EVm, h1, 'replicate');
    else
        % 模式0: 自适应DOPU滤波
        [EQmm] = vWinAvgFiltOpt_2_1(EQm, dopu_splitSpec_M, kRL, kRU);
        [EUmm] = vWinAvgFiltOpt_2_1(EUm, dopu_splitSpec_M, kRL, kRU);
        [EVmm] = vWinAvgFiltOpt_2_1(EVm, dopu_splitSpec_M, kRL, kRU);
    end

    % 调用核心算法，同时获得处理前后的完整结果
    [LA, PhR, cumLA, LA_raw, PhR_raw, cumLA_raw] = FreeSpace_PSOCT_3_DDG_rmBG_7(EQmm, EUmm, EVmm, Stru_E, test_seg_top, h1, h2, Avnum);
end
%% ====================================================================================
% 函数名: volFlatten (体积数据展平)
% 功能: 根据检测到的组织表面对3D/4D体积数据进行展平处理
% 输入参数:
%   LAs - 输入体积数据 [nZ×nX×nCh×nY]
%         nZ: 深度层数, nX: A-line数目, nCh: 通道数, nY: B-scan数目
%   topLines - 每个位置的表面深度 [nX×nY]，来自表面分割算法
% 输出参数:
%   LAs_flat - 展平后的体积数据 [nZ×nX×nCh×nY]
% 算法原理:
%   - 将每个A-line的数据向上移动，使组织表面对齐到统一深度
%   - 消除表面不平整造成的深度变化，便于后续分层分析
%   - 展平后的数据保持原始深度信息，但表面位置归一化
% 应用场景:
%   - 视网膜OCT分层分析
%   - 组织结构的定量比较
%   - 病理区域的深度校正
% ====================================================================================
function [LAs_flat] = volFlatten(LAs, topLines)
    % 获取输入数据的维度信息
    if ndims(LAs) == 4
        [nZ, nX, nCh, nY] = size(LAs);
        is4D = true;
    else
        [nZ, nX, nY] = size(LAs);
        nCh = 1;
        is4D = false;
    end

    % 初始化输出矩阵
    if is4D
        LAs_flat = zeros(nZ, nX, nCh, nY);
    else
        LAs_flat = zeros(nZ, nX, nY);
    end

    % 逐个A-line进行展平处理
    for i = 1:nX          % 遍历每条A-line
        for j = 1:nY      % 遍历每个B-scan
            % 计算从表面开始的有效数据长度
            % topLines(i,j)为该位置的表面深度索引
            effective_length = nZ - topLines(i,j) + 1;
            
            % 将表面以下的数据移动到顶部
            if effective_length > 0
                if is4D
                    % 4D情况：原始数据有多个通道
                    LAs_flat(1:effective_length, i, :, j) = ...
                        LAs(topLines(i,j):nZ, i, :, j);
                else
                    % 3D情况：单通道数据
                    LAs_flat(1:effective_length, i, j) = ...
                        LAs(topLines(i,j):nZ, i, j);
                end
            end
            
            % 注意：表面以上的区域在LAs_flat中保持为零
            % 这样可以清楚地区分组织内部和背景区域
        end
    end
end

%% ====================================================================================
% 函数名: quColoring (QU彩色编码)
% 功能: 将双折射的Q、U分量转换为HSV彩色编码图像
% 输入参数:
%   fLAs - 双折射数据矩阵 [nZ×nX×3×nY]，第3维为[?, Q, U]分量
%   cmapshift - 颜色映射旋转角度(可选参数)，用于调整色彩显示
% 输出参数:
%   hsvLA - HSV彩色编码图像 [nZ×nX×3×nY]，第3维为[H,S,V]分量
% 算法原理:
%   - 使用Q、U分量计算角度: θ = atan2(U, Q)
%   - 角度映射到色调(Hue): θ ∈ [-π, π] → H ∈ [0, 512]
%   - 使用HSV颜色空间表示双折射方向和强度
%   - 支持颜色映射的循环移位以优化显示效果
% 应用场景:
%   - 双折射方向的可视化
%   - 胶原纤维取向分析
%   - 病理组织的偏振特征显示
% 颜色含义:
%   - 色调(H): 双折射轴的方向角
%   - 饱和度(S): 通常设为最大值
%   - 亮度(V): 可用于表示双折射强度
% ====================================================================================
function [hsvLA] = quColoring(fLAs, cmapshift)
    % 创建HSV颜色映射表(512色)
    cusmapRg = hsv(512);

    % 备用颜色映射: cusmapRg = cusColormap(0);

    % 如果提供了颜色映射移位参数，则循环移位颜色表
    if nargin > 1
        cusmapRg = circshift(cusmapRg, cmapshift);
    end

    % 备用颜色映射: cusmapRg = [hsv(256); hsv(256)];

    % 定义角度范围[-π, π]对应的索引范围[1, 512]
    thetaRg = linspace(-pi, pi, 256*2);

    % 计算双折射轴角度: θ = atan2(U, Q)
    % fLAs(:,:,2) = U分量, fLAs(:,:,1) = Q分量
    thetas = atan2(fLAs(:,:,2), fLAs(:,:,1));


    thpf = polyfit(thetaRg,1:256*2,1);
    thetasInds = round(polyval(thpf,thetas));
    colorInds = cusmapRg(thetasInds,:);
    hsvLA = reshape(colorInds,[size(thetasInds),3]);
end

%% ====================================================================================
% 函数名: convert_dcm_to_tiff_local
% 功能: 快速提取DCM文件指定帧并转换为单帧TIFF文件
% 输入参数:
%   dcm_folder - DCM文件所在文件夹路径
%   target_frame - 目标帧号(从1开始)
%   output_prefix - 输出文件前缀名(可选)
% 输出: 在dcm_folder内创建tiff文件夹，保存单帧TIFF文件
% ====================================================================================
function convert_dcm_to_tiff_local(dcm_folder, target_frame, output_prefix)
    % 参数检查和默认值设置
    if nargin < 2
        error('至少需要提供DCM文件夹路径和目标帧号');
    end

    if nargin < 3
        output_prefix = '';
    end

    % 检查DCM文件夹是否存在
    if ~exist(dcm_folder, 'dir')
        error(['DCM文件夹不存在: ' dcm_folder]);
    end

    % 创建输出目录 - 直接在DCM文件夹内创建tiff子目录
    tiff_output_dir = fullfile(dcm_folder, 'tiff');

    if ~exist(tiff_output_dir, 'dir')
        mkdir(tiff_output_dir);
    end

    % 查找所有DCM文件
    dcm_files = dir(fullfile(dcm_folder, '*.dcm'));

    if isempty(dcm_files)
        warning(['在指定目录中未找到DCM文件: ' dcm_folder]);
        return;
    end

    fprintf('转换 %d 个DCM文件的第%d帧到TIFF...\n', length(dcm_files), target_frame);

    % 处理每个DCM文件
    success_count = 0;
    for i = 1:length(dcm_files)
        dcm_filename = dcm_files(i).name;
        dcm_filepath = fullfile(dcm_folder, dcm_filename);
        
        try
            % 读取DCM文件
            dcm_data = dicomread(dcm_filepath);
            
            % 确定帧的位置
            if ndims(dcm_data) == 4
                [~, ~, ~, total_frames] = size(dcm_data);
                frame_dim = 4;
            elseif ndims(dcm_data) == 3
                [~, ~, total_frames] = size(dcm_data);
                frame_dim = 3;
            else
                % 2D图像，直接保存
                frame_dim = 0;
                total_frames = 1;
            end
            
            % 提取目标帧
            if frame_dim == 0
                frame_data = dcm_data;
            else
                actual_frame = min(target_frame, total_frames);
                if frame_dim == 4
                    frame_data = dcm_data(:, :, :, actual_frame);
                else
                    frame_data = dcm_data(:, :, actual_frame);
                end
            end
            
            % 生成输出文件名
            [~, base_name, ~] = fileparts(dcm_filename);
            if isempty(output_prefix)
                tiff_filename = sprintf('%s_frame%d.tiff', base_name, target_frame);
            else
                tiff_filename = sprintf('%s_%s_frame%d.tiff', output_prefix, base_name, target_frame);
            end
            tiff_filepath = fullfile(tiff_output_dir, tiff_filename);
            
            % 保存为单帧TIFF文件
            imwrite(frame_data, tiff_filepath);
            success_count = success_count + 1;
            
        catch ME
            fprintf('处理文件 %s 出错: %s\n', dcm_filename, ME.message);
            continue;
        end
    end

    fprintf('完成! 成功转换 %d/%d 个文件到: %s\n', success_count, length(dcm_files), tiff_output_dir);
end
