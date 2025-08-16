%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% surface segmentation %%%%%%%%%%%%%%%%%%%%%%%%%
% ʹ��OAC���ݼ����Ʒ�����ȼ����ݶȣ�Ȼ�����һ����ֵ
% OAC: �����OAC����
% surf_threshold: �����������ֵ�ĵ�һ����ֵ

function test_seg_top=surf_seg(OAC,surf_threshold)

h1=[-0.2;-0.2;-0.2;-0.2;-0.2;1]; %�ϱ����ݶ�ģ��
OAC_G=(filter2(h1,OAC)>surf_threshold); %�ϱ���ģ���˲�
se = strel('disk', 2); % ѡ��һ���뾶Ϊ5��Բ�νṹԪ��
OAC_G = imdilate(OAC_G, se);
OAC_G = imerode(OAC_G, se);
OAC_G = double(bwareaopen(OAC_G, 20));
%figure,imshow(OAC_G,[0 1]);
for j=1:size(OAC,2)     
[~,locs] = findpeaks(OAC_G(:,j),'minpeakheight',0.2);% mph_ILM �趨��ֵ����С�߶�
        if isempty(locs) 
            locs(1)=Temp_locs;
        end
        ILM(j)=locs(1)-1; 
        Temp_locs=locs(1);
end
ILM = medfilt1(ILM,11);
test_seg_top = round(smooth(ILM, 15))+1;

end