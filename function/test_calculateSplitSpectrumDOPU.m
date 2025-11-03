%% ====================================================================================
% calculateSplitSpectrumDOPU 函数测试脚本
% 用于验证分裂谱DOPU计算函数的正确性
% ====================================================================================

% 添加函数路径
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
function_path = fullfile(parent_dir, 'function');
if exist(function_path, 'dir')
    addpath(function_path);
end

%% 1. 准备测试数据
fprintf('准备测试数据...\n');

% 模拟参数设置
SPL = 1024;          % FFT点数
nX = 10;             % A-scan数量 (缩小用于测试)
nr = 2;              % 重复次数
nWin = 9;            % 分裂谱窗口数

% 计算窗口参数
winL = floor(2*SPL/(nWin+1));           % 窗口长度
winG = tukeywin(winL, 0.25);            % 高斯窗口
windex = 1 : floor(winL/2) : SPL;       % 窗口起始位置

% 调整windex确保不超过范围
valid_indices = windex + winL - 1 <= SPL;
windex = windex(valid_indices);
nWin = length(windex);  % 更新实际窗口数

fprintf('窗口参数: winL=%d, nWin=%d\n', winL, nWin);
fprintf('windex = [%s]\n', num2str(windex));

% Z方向范围
czrg = 1:50;         % 感兴趣的深度范围 (缩小用于测试)
nZcrop = numel(czrg);

% 生成模拟OCT信号 (两个偏振通道)
% Bd1, Bd2: [信号长度 × X像素 × 重复次数]
Bd1 = randn(SPL, nX, nr) + 1i * randn(SPL, nX, nr);
Bd2 = randn(SPL, nX, nr) + 1i * randn(SPL, nX, nr);

%% 2. 设置参数结构体
params.dopu.do_ssdopu = true;  % 启用分裂谱计算

%% 3. 调用函数进行测试
fprintf('调用 calculateSplitSpectrumDOPU 函数...\n');

try
    tic;
    [dopu_splitSpectrum, dopu_ss] = calculateSplitSpectrumDOPU(...
        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);
    elapsed_time = toc;

    fprintf('函数调用成功，耗时: %.3f 秒\n', elapsed_time);

    %% 4. 验证结果
    fprintf('验证结果...\n');

    % 检查输出尺寸
    expected_size_dopu = [nZcrop, nX];
    expected_size_dopu_ss = [nZcrop, nX, nr];

    if size(dopu_splitSpectrum) == expected_size_dopu
        fprintf('✓ dopu_splitSpectrum 尺寸正确: [%d, %d]\n', size(dopu_splitSpectrum, 1), size(dopu_splitSpectrum, 2));
    else
        fprintf('✗ dopu_splitSpectrum 尺寸错误，期望: [%d, %d]，实际: [%d, %d]\n', ...
            expected_size_dopu(1), expected_size_dopu(2), size(dopu_splitSpectrum, 1), size(dopu_splitSpectrum, 2));
    end

    if size(dopu_ss) == expected_size_dopu_ss
        fprintf('✓ dopu_ss 尺寸正确: [%d, %d, %d]\n', size(dopu_ss, 1), size(dopu_ss, 2), size(dopu_ss, 3));
    else
        fprintf('✗ dopu_ss 尺寸错误，期望: [%d, %d, %d]，实际: [%d, %d, %d]\n', ...
            expected_size_dopu_ss(1), expected_size_dopu_ss(2), expected_size_dopu_ss(3), ...
            size(dopu_ss, 1), size(dopu_ss, 2), size(dopu_ss, 3));
    end

    % 检查数值范围 (DOPU应该在0-1之间)
    if all(dopu_splitSpectrum(:) >= 0) && all(dopu_splitSpectrum(:) <= 1)
        fprintf('✓ dopu_splitSpectrum 数值范围正确: [0, 1]\n');
    else
        fprintf('✗ dopu_splitSpectrum 数值范围错误\n');
        fprintf('  最小值: %.3f, 最大值: %.3f\n', min(dopu_splitSpectrum(:)), max(dopu_splitSpectrum(:)));
    end

    if all(dopu_ss(:) >= 0) && all(dopu_ss(:) <= 1)
        fprintf('✓ dopu_ss 数值范围正确: [0, 1]\n');
    else
        fprintf('✗ dopu_ss 数值范围错误\n');
        fprintf('  最小值: %.3f, 最大值: %.3f\n', min(dopu_ss(:)), max(dopu_ss(:)));
    end

    % 显示统计信息
    fprintf('统计信息:\n');
    fprintf('  dopu_splitSpectrum - 均值: %.3f, 标准差: %.3f\n', ...
        mean(dopu_splitSpectrum(:)), std(dopu_splitSpectrum(:)));
    fprintf('  dopu_ss - 均值: %.3f, 标准差: %.3f\n', ...
        mean(dopu_ss(:)), std(dopu_ss(:)));

    fprintf('测试完成！\n');

catch ME
    fprintf('测试失败: %s\n', ME.message);
    fprintf('错误详情:\n');
    disp(ME);
end

%% 5. 测试不进行分裂谱计算的情况
fprintf('\n测试不进行分裂谱计算的情况...\n');
params.dopu.do_ssdopu = false;

try
    [dopu_disabled, dopu_ss_disabled] = calculateSplitSpectrumDOPU(...
        Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);

    if all(dopu_disabled(:) == 1)
        fprintf('✓ 禁用分裂谱计算时返回正确默认值 (全1矩阵)\n');
    else
        fprintf('✗ 禁用分裂谱计算时返回值错误\n');
    end

catch ME
    fprintf('禁用测试失败: %s\n', ME.message);
end

fprintf('\n所有测试完成！\n');