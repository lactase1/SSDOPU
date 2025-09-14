% GPU利用率优化测试脚本
% 测试不同GPU优化策略的效果

clear all; clc;

fprintf('==================================================\n');
fprintf('GPU利用率优化验证测试\n');
fprintf('==================================================\n');

% 检测GPU
try
    gpu_info = gpuDevice();
    fprintf('GPU: %s\n', gpu_info.Name);
    fprintf('内存: %.2f GB / %.2f GB\n', gpu_info.AvailableMemory/1024^3, gpu_info.TotalMemory/1024^3);
    fprintf('多处理器: %d\n', gpu_info.MultiprocessorCount);
catch
    error('未检测到GPU，无法进行测试');
end

% 模拟PSOCT数据尺寸
fprintf('\n==================================================\n');
fprintf('创建测试数据...\n');
fprintf('==================================================\n');

% 模拟实际PSOCT数据大小
Blength = 4096;  % 光谱长度
nX = 512;        % A-line数量
nr = 7;          % 重复次数
nWin = 9;        % split spectrum窗口数
winL = 2*Blength/(nWin+1);
SPL = 4096;
nZcrop = 300;

fprintf('数据尺寸: %d x %d x %d\n', Blength, nX, nr);
fprintf('内存需求: 约%.2f GB\n', Blength*nX*nr*8*2/1024^3); % 复数双精度

%% 测试1: 原始方法（低GPU利用率）
fprintf('\n==================================================\n');
fprintf('测试1: 原始方法 (模拟低GPU利用率)\n');
fprintf('==================================================\n');

% 创建测试数据
test_data1 = randn(Blength, nX, nr) + 1i*randn(Blength, nX, nr);
test_data2 = randn(Blength, nX, nr) + 1i*randn(Blength, nX, nr);
winG = tukeywin(winL, 0.25);
windex = 1:winL/2:Blength;

tic;
% 原始方法：频繁的小块GPU传输
for simulation = 1:3  % 模拟多次调用
    % 每次都重新传输数据（低效）
    gpu_data1 = gpuArray(test_data1);
    gpu_data2 = gpuArray(test_data2);
    winG_gpu = gpuArray(winG);
    
    % 嵌套循环（低GPU利用率）
    Bimg1 = zeros(SPL, nX, nr, nWin, 'gpuArray');
    Bimg2 = zeros(SPL, nX, nr, nWin, 'gpuArray');
    
    for iR = 1:nr
        for iL = 1:nWin
            % 小块操作
            data_fragment1 = gpu_data1(windex(iL):windex(iL)+winL-1, :, iR) .* winG_gpu;
            data_fragment2 = gpu_data2(windex(iL):windex(iL)+winL-1, :, iR) .* winG_gpu;
            Bimg1(:,:,iR,iL) = fft(data_fragment1, SPL, 1);
            Bimg2(:,:,iR,iL) = fft(data_fragment2, SPL, 1);
        end
    end
    
    % 频繁的gather操作
    result_low = gather(mean(abs(Bimg1).^2 + abs(Bimg2).^2, [3,4]));
    clear gpu_data1 gpu_data2 winG_gpu Bimg1 Bimg2;
end
time_low = toc;
fprintf('原始方法时间: %.3f 秒\n', time_low);

%% 测试2: 优化方法（高GPU利用率）
fprintf('\n==================================================\n');
fprintf('测试2: 优化方法 (高GPU利用率)\n');
fprintf('==================================================\n');

tic;
% 优化方法：减少数据传输，批量化操作
% 一次性传输静态数据
winG_gpu_opt = gpuArray(winG);
winG_whole_gpu = gpuArray(tukeywin(Blength, 0.25));

for simulation = 1:3
    % 批量传输数据
    gpu_data1_batch = gpuArray(test_data1);
    gpu_data2_batch = gpuArray(test_data2);
    
    % 预分配大块GPU内存
    all_results = zeros(SPL, nX, nr, nWin, 'gpuArray');
    
    % 批量化FFT计算
    for iR = 1:nr
        % 一次处理所有窗口
        windowed_batch1 = zeros(winL, nX, nWin, 'gpuArray');
        windowed_batch2 = zeros(winL, nX, nWin, 'gpuArray');
        
        for iL = 1:nWin
            windowed_batch1(:,:,iL) = gpu_data1_batch(windex(iL):windex(iL)+winL-1,:,iR) .* winG_gpu_opt;
            windowed_batch2(:,:,iL) = gpu_data2_batch(windex(iL):windex(iL)+winL-1,:,iR) .* winG_gpu_opt;
        end
        
        % 批量FFT
        fft_batch1 = fft(windowed_batch1, SPL, 1);
        fft_batch2 = fft(windowed_batch2, SPL, 1);
        
        % 批量化Stokes参数计算
        S0_batch = abs(fft_batch1).^2 + abs(fft_batch2).^2;
        all_results(:,:,iR,:) = S0_batch;
    end
    
    % 减少gather次数
    result_high = gather(mean(all_results, [3,4]));
    clear gpu_data1_batch gpu_data2_batch all_results;
end
time_high = toc;
fprintf('优化方法时间: %.3f 秒\n', time_high);

%% 结果对比
fprintf('\n==================================================\n');
fprintf('GPU利用率优化效果\n');
fprintf('==================================================\n');

speedup = time_low / time_high;
efficiency_improvement = (1 - time_high/time_low) * 100;

fprintf('原始方法: %.3f 秒\n', time_low);
fprintf('优化方法: %.3f 秒\n', time_high);
fprintf('加速比: %.2fx\n', speedup);
fprintf('效率提升: %.1f%%\n', efficiency_improvement);

% GPU利用率提升分析
fprintf('\n优化策略效果分析:\n');
if speedup > 2
    fprintf('✅ 显著提升! GPU利用率优化效果很好\n');
    fprintf('   - 减少了数据传输开销\n');
    fprintf('   - 提高了GPU并行度\n');
    fprintf('   - 批量化操作更高效\n');
elseif speedup > 1.3
    fprintf('✅ 有明显提升! GPU利用率得到改善\n');
    fprintf('   - 建议进一步优化内存访问模式\n');
else
    fprintf('⚠️  提升有限，可能需要:\n');
    fprintf('   - 检查数据传输瓶颈\n');
    fprintf('   - 进一步增加批量大小\n');
    fprintf('   - 优化内存访问模式\n');
end

%% GPU内存使用效率测试
fprintf('\n==================================================\n');
fprintf('GPU内存使用效率\n');
fprintf('==================================================\n');

% 测试大块内存分配 vs 小块分配
fprintf('测试内存分配策略...\n');

% 小块分配（低效）
tic;
small_arrays = cell(100, 1);
for i = 1:100
    small_arrays{i} = zeros(512, 512, 'gpuArray');
end
clear small_arrays;
time_small = toc;

% 大块分配（高效）
tic;
large_array = zeros(512, 512, 100, 'gpuArray');
clear large_array;
time_large = toc;

fprintf('小块分配: %.4f 秒\n', time_small);
fprintf('大块分配: %.4f 秒\n', time_large);
fprintf('内存效率提升: %.1fx\n', time_small/time_large);

%% 清理GPU
reset(gpuDevice);
fprintf('\n🧹 GPU内存已清理\n');

fprintf('\n==================================================\n');
fprintf('测试完成!\n');
fprintf('==================================================\n');

if speedup > 1.5
    fprintf('🎉 GPU利用率优化成功!\n');
    fprintf('   现在应该能看到更高的GPU利用率了。\n');
else
    fprintf('🔧 需要进一步优化GPU利用率\n');
    fprintf('   建议检查数据大小和并行策略。\n');
end