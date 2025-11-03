% 自动读取 tmp 文件夹下所有 mat 文件并比较格式
mat_dir = 'D:/1-Liu Jian/yongxin.wang/tmp/';
files = dir(fullfile(mat_dir, '*.mat'));

if numel(files) < 2
    error('mat 文件少于两个，无法比较！');
end

for i = 1:numel(files)
    fprintf('第 %d 个文件: %s\n', i, files(i).name);
    info = whos('-file', fullfile(mat_dir, files(i).name));
    disp(info);
end

% 只比较前两个文件的 topLines 格式
file1 = fullfile(mat_dir, files(1).name);
file2 = fullfile(mat_dir, files(2).name);
data1 = load(file1);
data2 = load(file2);

if isfield(data1, 'topLines') && isfield(data2, 'topLines')
    disp('topLines 尺寸对比：');
    disp(['file1: ', mat2str(size(data1.topLines)), ', class: ', class(data1.topLines)]);
    disp(['file2: ', mat2str(size(data2.topLines)), ', class: ', class(data2.topLines)]);
else
    disp('至少有一个文件不包含 topLines 变量');
end
