%% ====================================================================================
% PS-OCT 分裂谱DOPU计算模块 - 使用指南
% 版本: 1.0
% 更新日期: 2024年10月18日
% ====================================================================================

%% 概述
% 本模块将PS-OCT处理流程中的分裂谱DOPU计算部分抽离为独立函数，
% 便于代码复用、测试和维护。

%% 文件结构
% function/
%   ├── calculateSplitSpectrumDOPU.m      % 主函数：分裂谱DOPU计算
%   ├── cumulativeQUV.m                   % 辅助函数：Stokes参数计算
%   ├── test_calculateSplitSpectrumDOPU.m % 测试脚本
%   └── calculateSplitSpectrumDOPU_README.m % 详细文档

%% 快速开始

%% 1. 基本用法
% 在PS-OCT处理脚本中替换原有的分裂谱代码：

% 原来的代码 (已注释):
% if ~params.dopu.do_ssdopu
%     dopu_ss = 1;
% else
%     % ... 复杂的分裂谱计算代码 ...
%     dopu_splitSpectrum(:,:,iY) = mean(dopu_ss,3);
% end

% 替换为:
[dopu_splitSpectrum(:,:,iY), dopu_ss] = calculateSplitSpectrumDOPU(...
    Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);

%% 2. 参数准备
% 确保以下参数已正确设置:

% 基本参数
SPL = 1024;           % FFT点数
nX = 512;             % X方向像素数
nr = 4;               % 重复次数
nWin = 9;             % 分裂谱窗口数

% 窗口参数计算
winL = 2*SPL/(nWin+1);                    % 窗口长度
winG = tukeywin(winL, 0.25);              % 高斯窗口
windex = 1 : winL/2 : SPL;                % 窗口起始位置

% Z方向裁剪
czrg = 1:320;          % 感兴趣的深度范围

% 参数结构体
params.dopu.do_ssdopu = true;  % 启用分裂谱计算

%% 3. 完整调用示例
% % 准备输入数据 (Bd1, Bd2 应已进行背景减除和色散校正)
% [dopu_result, dopu_all_reps] = calculateSplitSpectrumDOPU(...
%     Bd1, Bd2, params, SPL, nX, nr, nWin, windex, winL, winG, czrg);
%
% % dopu_result: [nZcrop × nX] 最终DOPU结果
% % dopu_all_reps: [nZcrop × nX × nr] 各重复的DOPU结果

%% 4. 测试函数
% 运行测试验证功能:
run('function/test_calculateSplitSpectrumDOPU.m');

%% 算法说明

%% 分裂谱DOPU计算流程:
% 1. 将OCT频谱分割成多个重叠子频带 (默认9个)
% 2. 对每个子频带应用高斯窗口并进行FFT
% 3. 计算每个子频带的Stokes参数 (S0, S1, S2, S3)
% 4. 对所有子频带进行平均
% 5. 计算DOPU = sqrt(<S1/S0>² + <S2/S0>² + <S3/S0>²)

%% 输出说明
% - dopu_splitSpectrum: 最终DOPU图像 [Z × X]
% - dopu_ss: 各重复次数的DOPU [Z × X × 重复次数]
% - DOPU取值范围: [0, 1]，1表示完全偏振，0表示非偏振

%% 性能优化
% - 支持并行处理 (parfor循环)
% - 内存预分配以提高性能
% - 自动参数验证和错误处理

%% 依赖关系
% - MATLAB Signal Processing Toolbox (tukeywin, fft)
% - cumulativeQUV.m (Stokes参数计算)

%% 版本历史
% v1.0 (2024-10-18)
% - 初始版本，从主处理脚本中抽离
% - 支持标准分裂谱DOPU计算
% - 包含完整的参数验证和错误处理
% - 提供详细的文档和测试

%% 故障排除

%% 常见问题:
% 1. "未定义函数 cumulativeQUV"
%    解决: 确保 function/ 目录在MATLAB路径中
%
% 2. 内存不足错误
%    解决: 减小 nWin 或使用更小的测试数据集
%
% 3. 尺寸不匹配错误
%    解决: 检查输入参数的尺寸是否与预期一致

%% 技术支持
% 如有问题，请检查:
% 1. 函数文档: calculateSplitSpectrumDOPU_README.m
% 2. 测试脚本: test_calculateSplitSpectrumDOPU.m
% 3. 原始代码: rPSOCT_06_splitSpec_3_3D_without_mat_wyx.m

%% ====================================================================================