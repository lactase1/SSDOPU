function print_progress(current, total, prefix, bar_length, style)
% print_progress - 在MATLAB命令窗口显示纯文本进度条
% 
% 输入参数:
%   current    - 当前进度值 (1 到 total)
%   total      - 总数
%   prefix     - 进度条前缀文字 (可选，默认为 'Progress')
%   bar_length - 进度条长度 (可选，默认为 30)
%   style      - 进度条样式 (可选，默认为 'gradient')
%                'gradient' - 渐变符号样式 (#→*→=)
%                'blocks'   - 方块样式 (█)
%                'arrows'   - 箭头样式 (>)
%                'simple'   - 简单样式 (=)
%
% 示例:
%   for i = 1:100
%       print_progress(i, 100, '处理中', 30, 'gradient');
%       pause(0.05);
%   end

    % 参数默认值
    if nargin < 3 || isempty(prefix)
        prefix = 'Progress';
    end
    if nargin < 4 || isempty(bar_length)
        bar_length = 30;
    end
    if nargin < 5 || isempty(style)
        style = 'gradient';
    end
    
    % 计算进度百分比
    percent = current / total;
    filled_length = round(bar_length * percent);
    
    % 根据样式构建进度条（纯文本版）
    switch lower(style)
        case 'gradient'
            bar = build_gradient_bar_text(filled_length, bar_length, percent);
        case 'blocks'
            bar = build_blocks_bar_text(filled_length, bar_length);
        case 'arrows'
            bar = build_arrows_bar_text(filled_length, bar_length);
        case 'simple'
            bar = build_simple_bar_text(filled_length, bar_length);
        otherwise
            bar = build_gradient_bar_text(filled_length, bar_length, percent);
    end
    
    % 打印进度条（使用 \r 回到行首实现原地更新）
    if current < total
        fprintf('\r%s [%s] %5.1f%% (%d/%d)', ...
                prefix, bar, percent*100, current, total);
    else
        % 完成时换行并显示特殊完成标记
        fprintf('\r%s [%s] 100.0%% (%d/%d) ✓\n', ...
                prefix, bar, current, total);
    end
    
    % 强制刷新输出
    drawnow;
end

%% ============ 样式构建函数 ============

%% ============ 辅助函数（纯文本版） ============

function bar = build_gradient_bar_text(filled, total, percent)
    % 渐变符号进度条（统一使用=符号）
    if filled > 0
        filled_part = repmat('=', 1, filled);
    else
        filled_part = '';
    end
    
    empty_part = repmat('.', 1, total - filled);
    bar = [filled_part empty_part];
end

function bar = build_blocks_bar_text(filled, total)
    % 方块样式进度条（纯文本）
    filled_part = repmat('█', 1, filled);
    empty_part = repmat('░', 1, total - filled);
    bar = [filled_part empty_part];
end

function bar = build_arrows_bar_text(filled, total)
    % 箭头样式进度条（纯文本）
    if filled > 0
        filled_part = repmat('>', 1, filled);
    else
        filled_part = '';
    end
    empty_part = repmat('-', 1, total - filled);
    bar = [filled_part empty_part];
end

function bar = build_simple_bar_text(filled, total)
    % 简单样式（纯文本）
    filled_part = repmat('=', 1, filled);
    if filled < total
        arrow = '>';
    else
        arrow = '';
    end
    empty_part = repmat(' ', 1, max(0, total - filled - 1));
    bar = [filled_part arrow empty_part];
end

function str = format_time(seconds)
    % 将秒数格式化为可读的时间字符串
    if isnan(seconds) || isinf(seconds)
        str = 'N/A';
    elseif seconds < 60
        str = sprintf('%.0f秒', seconds);
    elseif seconds < 3600
        str = sprintf('%.0f分%.0f秒', floor(seconds/60), mod(seconds, 60));
    else
        str = sprintf('%.0f时%.0f分', floor(seconds/3600), mod(floor(seconds/60), 60));
    end
end
