%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% surface segmentation %%%%%%%%%%%%%%%%%%%%%%%%%
% 使用OAC数据检测样品，首先计算梯度，然后检测第一个峰值
% OAC: 输入的OAC数据
% surf_threshold: 检测大于这个阈值的第一个峰值

function test_seg_top=surf_seg(OAC,surf_threshold)

h1=[-0.2;-0.2;-0.2;-0.2;-0.2;1]; %上表面梯度模版
OAC_G=(filter2(h1,OAC)>surf_threshold); %上表面模板滤波
se = strel('disk', 2); % 选择一个半径为5的圆形结构元素
OAC_G = imdilate(OAC_G, se);
OAC_G = imerode(OAC_G, se);
OAC_G = double(bwareaopen(OAC_G, 20));
%figure,imshow(OAC_G,[0 1]);
for j=1:size(OAC,2)     
[~,locs] = findpeaks(OAC_G(:,j),'minpeakheight',0.2);% mph_ILM 设定峰值的最小高度
        if isempty(locs) 
            locs(1)=Temp_locs;
        end
        ILM(j)=locs(1)-1; 
        Temp_locs=locs(1);
end
ILM = medfilt1(ILM,11);
test_seg_top = round(smooth(ILM, 15))+1;

end