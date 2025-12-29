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
