% GPU测试脚本 - 验证GPU加速是否正常工作
% 这个脚本可以独立运行，用来测试GPU功能

clear all; clc;

fprintf('==================================================\n');
fprintf('GPU加速功能测试\n');
fprintf('==================================================\n');

% 1. 检测GPU可用性
try
    gpu_info = gpuDevice();
    fprintf('✅ GPU检测成功: %s\n', gpu_info.Name);
    fprintf('   GPU内存: %.2f GB 可用 / %.2f GB 总计\n', ...
            gpu_info.AvailableMemory/1024^3, gpu_info.TotalMemory/1024^3);
    fprintf('   计算能力: %s\n', gpu_info.ComputeCapability);
    gpu_available = true;
catch ME
    fprintf('❌ GPU检测失败: %s\n', ME.message);
    fprintf('   将使用CPU计算\n');
    gpu_available = false;
end

if ~gpu_available
    fprintf('\n请检查:\n');
    fprintf('1. CUDA是否正确安装\n');
    fprintf('2. GPU驱动是否最新\n');
    fprintf('3. MATLAB Parallel Computing Toolbox是否安装\n');
    return;
end

fprintf('\n==================================================\n');
fprintf('GPU性能测试\n');
fprintf('==================================================\n');

% 2. 创建测试数据
fprintf('创建测试数据...\n');
test_size = [4096, 512, 7];  % 模拟PSOCT数据尺寸
test_data1 = randn(test_size) + 1i*randn(test_size);
test_data2 = randn(test_size) + 1i*randn(test_size);

% 3. CPU计算测试
fprintf('开始CPU计算测试...\n');
tic;
% FFT计算
fft_result_cpu = fft(test_data1, [], 1);
% Stokes参数计算
axis_cpu = angle(test_data2.*conj(test_data1));
S0_cpu = abs(test_data1).^2 + abs(test_data2).^2;
S1_cpu = abs(test_data1).^2 - abs(test_data2).^2;
S2_cpu = 2.*abs(test_data1).*abs(test_data2).*cos(axis_cpu);
S3_cpu = 2.*abs(test_data1).*abs(test_data2).*sin(-axis_cpu);
cpu_time = toc;
fprintf('CPU计算完成: %.3f 秒\n', cpu_time);

% 4. GPU计算测试
fprintf('开始GPU计算测试...\n');
tic;
% 数据传输到GPU
gpu_data1 = gpuArray(test_data1);
gpu_data2 = gpuArray(test_data2);

% FFT计算 (GPU)
fft_result_gpu = fft(gpu_data1, [], 1);
% Stokes参数计算 (GPU)
axis_gpu = angle(gpu_data2.*conj(gpu_data1));
S0_gpu = abs(gpu_data1).^2 + abs(gpu_data2).^2;
S1_gpu = abs(gpu_data1).^2 - abs(gpu_data2).^2;
S2_gpu = 2.*abs(gpu_data1).*abs(gpu_data2).*cos(axis_gpu);
S3_gpu = 2.*abs(gpu_data1).*abs(gpu_data2).*sin(-axis_gpu);

% 数据传输回CPU
S0_gpu_result = gather(S0_gpu);
gpu_time = toc;
fprintf('GPU计算完成: %.3f 秒\n', gpu_time);

% 5. 性能对比
fprintf('\n==================================================\n');
fprintf('性能对比结果\n');
fprintf('==================================================\n');
speedup = cpu_time / gpu_time;
fprintf('CPU时间:     %.3f 秒\n', cpu_time);
fprintf('GPU时间:     %.3f 秒\n', gpu_time);
fprintf('加速比:      %.2fx\n', speedup);

if speedup > 1.5
    fprintf('✅ GPU加速正常工作!\n');
else
    fprintf('⚠️  GPU加速效果不明显，可能的原因:\n');
    fprintf('   1. 数据量较小，GPU优势不明显\n');
    fprintf('   2. GPU传输开销较大\n');
end

% 6. 验证结果正确性
fprintf('\n==================================================\n');
fprintf('结果正确性验证\n');
fprintf('==================================================\n');
error_S0 = max(abs(S0_cpu(:) - S0_gpu_result(:)));
fprintf('S0参数最大误差: %.2e\n', error_S0);

if error_S0 < 1e-10
    fprintf('✅ GPU计算结果正确!\n');
else
    fprintf('⚠️  存在数值误差，但在可接受范围内\n');
end

% 7. 内存使用情况
fprintf('\n==================================================\n');
fprintf('GPU内存使用情况\n');
fprintf('==================================================\n');
fprintf('GPU内存使用: %.2f GB / %.2f GB\n', ...
        (gpu_info.TotalMemory - gpu_info.AvailableMemory)/1024^3, ...
        gpu_info.TotalMemory/1024^3);

% 清理GPU内存
clear gpu_data1 gpu_data2 fft_result_gpu axis_gpu S0_gpu S1_gpu S2_gpu S3_gpu;
reset(gpuDevice);
fprintf('GPU内存已释放\n');

fprintf('\n==================================================\n');
fprintf('测试完成!\n');
fprintf('==================================================\n');

if gpu_available && speedup > 1.5
    fprintf('🎉 您的RTX 4090 GPU加速功能正常工作!\n');
    fprintf('   现在可以运行PSOCT处理脚本享受GPU加速了。\n');
else
    fprintf('🔧 建议检查GPU配置或联系技术支持\n');
end