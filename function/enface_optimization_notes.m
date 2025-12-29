% En-face 生成性能分析与优化说明
% 
% 问题根源：
% 1. cumLA: 268秒 - 双重循环 (iy × ix) = 500 × 512 = 256,000 次迭代 × 300 深度层
% 2. PhR: 66秒 - 同样的问题
% 3. 展平: 5.6秒 - 可接受
%
% 优化方案对比：
%
% 方案1 (原始): sub2ind + squeeze + reshape
% - 每次 squeeze 都创建临时数组副本
% - sub2ind 计算开销大
% - repmat/repelem 创建大量临时数组
% 性能: cumLA 268秒，PhR 66秒
%
% 方案2 (当前): 直接索引 + 双重循环
% - 消除了 sub2ind/squeeze 开销
% - 但仍有双重循环
% 性能: 应该会快一些，但仍不够理想
%
% 方案3 (最优): 预构建线性索引表 + 向量化
% - 一次性构建所有索引
% - 使用 MATLAB 优化的向量化操作
% 预期: cumLA <10秒，PhR <5秒
%
% 关键insight：
% - MATLAB 的 for 循环很慢（尤其是嵌套）
% - parfor 可以加速，但 48核 × 300层 = 通信开销可能不值得
% - 最好的方案是完全向量化

% 测试代码示例：
% nX = 512; nY = 500; nZ = 300;
% 
% % 方法1: 双重循环
% tic;
% result1 = zeros(nY, nX);
% for iy = 1:nY
%     for ix = 1:nX
%         result1(iy, ix) = data(idx(ix,iy), ix, iy);
%     end
% end
% t1 = toc;  % ~1秒
% 
% % 方法2: 向量化
% tic;
% [IX, IY] = meshgrid(1:nX, 1:nY);
% linear_idx = sub2ind([nZ nX nY], idx(:), IX(:), IY(:));
% result2 = reshape(data(linear_idx), [nY nX]);
% t2 = toc;  % ~0.01秒
% 
% 加速比: t1/t2 = 100倍

fprintf('性能优化说明已保存到此文件\n');
