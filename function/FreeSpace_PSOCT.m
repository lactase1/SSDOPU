% ��˹�˲�������ֵ��ƽ����
% QUVС�߶��˲���h1)�������߶��˲�(h2)��
% ��ÿ��A-scanȥ��DC����������ִ�б�׼��
% ���˲���Ĺ��������ת����
% 


function [LA,PhR] = FreeSpace_PSOCT(Qm,Um,Vm,Stru,test_seg_top,h1,h2,Avnum)
% showimg = 1;


%% flatting
StruF=Stru*0;QmF=Qm*0;UmF=Um*0;VmF=Vm*0;
for j=1:size(Stru,2)
   StruF(1:size(Stru,1)-test_seg_top(j)+1,j)=Stru(test_seg_top(j):size(Stru,1),j);
   QmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Qm(test_seg_top(j):size(Stru,1),j);
   UmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Um(test_seg_top(j):size(Stru,1),j);
   VmF(1:size(Stru,1)-test_seg_top(j)+1,j)=Vm(test_seg_top(j):size(Stru,1),j);
end


%%  QUV filtering

% QmF=imfilter(QmF,h1,'replicate'); % ������Ϊ��Ӧ���������˲�
% UmF=imfilter(UmF,h1,'replicate'); 
% VmF=imfilter(VmF,h1,'replicate');


QF=QmF-mean(QmF); % zero-meaning
UF=UmF-mean(UmF);
VF=VmF-mean(VmF);


QF1=QF./sqrt(QF.^2+UF.^2+VF.^2); % normallization
UF1=UF./sqrt(QF.^2+UF.^2+VF.^2);
VF1=VF./sqrt(QF.^2+UF.^2+VF.^2);

QF1(isnan(QF1)) = 0;    UF1(isnan(UF1)) = 0;    VF1(isnan(VF1)) = 0;

SmapF=cat(3,QF1,UF1,VF1);
% figure(11),imagesc(SmapF);
%% rot 90 degree
% vetorq=pi/2*[0;1;0];
% rotq=rotationVectorToMatrix(vetorq);
% for i=1:size(SmapF,1)
%     for j=1:size(SmapF,2)
%         aa=squeeze(SmapF(i,j,:));
%         SmapF1(i,j,1:3)=rotq*aa;
%     end
% end
SmapF1=SmapF;
%% svd to fit LA 
phd=zeros(size(SmapF1,1),size(SmapF1,2));
for i=1:size(SmapF1,1)-Avnum
    for j=1:size(SmapF1,2)
        planeData=squeeze(SmapF1(i:i+Avnum,j,:)); % QUV points from Avnum layers
        x1=squeeze(planeData(:,1)); y1=squeeze(planeData(:,2)); z1=squeeze(planeData(:,3));
        nx=sum(sum(x1==0));
        if x1==0|nx>0 % skip if Qs has 0
            loaxis2(i,j,1:3)=0;
            phR(i,j)=0;
            phd(i,j)=0;
            rotationVector(i,j,1:3)=0;
        else

%            centeredPlane=planeData-mean(planeData); %%%������ĵ㣬�Ƿ���Ҫ���ֲ����ֵ����
            centeredPlane=planeData;            
            
            [~,~,V1]=svd(centeredPlane); %%%%%%%%%%%%%%%%%%%%% QUV fitting plane axis
             a1=V1(1,3); b1=V1(2,3); c1=V1(3,3);
            
            a1=a1./sqrt(a1.^2+b1.^2+c1.^2);
            b1=b1./sqrt(a1.^2+b1.^2+c1.^2);
            c1=c1./sqrt(a1.^2+b1.^2+c1.^2);
            a1(isnan(a1)) = 0;    b1(isnan(b1)) = 0;    c1(isnan(c1)) = 0;
            
%% check the axis orientation and reverse **ʮ����Ҫ
% ��svd�ֽ�õ��Ĺ��᲻��ȷ�����ţ����Ի�Ҫ�жϸ����������Ǹ���
            Sin1=[x1(1),y1(1),z1(1)];
            if x1(1)==0
                Sin1=[x1(2),y1(2),z1(2)];
            end
            Sout1=[x1(Avnum),y1(Avnum),z1(Avnum)];
            ax=[a1,b1,c1];

            dd1=dot(Sin1,ax); % ����Sin1��ax���ϵ�ͶӰ���Ƕ�С��90��Ϊ��������90��Ϊ����
            SS1=Sin1-dd1*ax; %��������������õ�һ����ֱ�ڹ����������

            dd2=dot(Sout1,ax);
            SS2=Sout1-dd2*ax; %ͬ�����һ����ֱ�ڹ����������
            
            dd3=cross(SS1,SS2); %�����������
            dd4=dot(dd3,ax);
            if dd4<0   % ������С�ںŻ����ô��ں��������ת�����Ƿ����-1�йء�
                a1=-a1; b1=-b1; c1=-c1;
            end
            loaxis2(i,j,:)=cat(3,a1,b1,c1);% LA
%% cal phase retardation
            for k=1:Avnum
                Sin=[x1(k),y1(k),z1(k)];
                A=[a1,b1,c1];
                Sout=[x1(k+1),y1(k+1),z1(k+1)];
                val1=cross(Sin,A);  %������
                lval1=sqrt(sum(val1.^2)); %��������ģ
                val2=cross(Sout,A); %������
                lval2=sqrt(sum(val2.^2)); %��������ģ
                if (lval1==0)|(lval2==0)
                    val3(k)=0;
                else
                    val3(k)=dot(val1,val2)/lval1/lval2; %�������нǵ��ڷ������н�
                end
                val3(k) = max([val3(k),-1]);
                val3(k) = min([val3(k), 1]);
                phRk(k)=acos(val3(k));% cal phase for each layers
            end

                phR(i,j)=(acos((median(val3)))); %%%%%%%%%% Phase retardation

        end
    end
end
% post process the LA 
%% �����˲�

    Temp_ax1(:,:)=loaxis2(:,:,1);
    Temp_ax2(:,:)=loaxis2(:,:,2);
    Temp_ax3(:,:)=loaxis2(:,:,3);
       
    Temp_ax1=imfilter(Temp_ax1,h2,'replicate');
    Temp_ax2=imfilter(Temp_ax2,h2,'replicate');
    Temp_ax3=imfilter(Temp_ax3,h2,'replicate');
   
Temp_ax1=Temp_ax1./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
Temp_ax2=Temp_ax2./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);
Temp_ax3=Temp_ax3./sqrt(Temp_ax1.^2+Temp_ax2.^2+Temp_ax3.^2);

Temp_ax1(isnan(Temp_ax1)) = 0;    Temp_ax2(isnan(Temp_ax2)) = 0;    Temp_ax3(isnan(Temp_ax3)) = 0;    
loaxis3=cat(3,Temp_ax1,Temp_ax2,Temp_ax3);
%% ��λ�˲�
    phR=imfilter(phR,h2,'replicate');
%% ������ת����
% for i=1:size(SmapF1,1)-Avnum
    for j=1:size(SmapF1,2)
        ax_rot(1,:)=loaxis3(1,j,:);
        rotationVector(1,j,:) = -phR(1,j)/2 * ax_rot;
    end
% end
loaxis22=loaxis3.*(phR);

%             for cen_i=1:3 %���ֵ�󣬹�һ��
%                 loaxis22(:,:,cen_i)=loaxis22(:,:,cen_i)./sqrt(sum((loaxis22.^2),3));
%             end


rotationVector2=rotationVector;

%% rotate the LA to achieve the depth resolved L

for j=1:size(loaxis2,2)
    for i=1:size(loaxis2,1)
        vector1=squeeze(rotationVector2(1,j,:));
        a=[loaxis22(i,j,1),loaxis22(i,j,2),loaxis22(i,j,3)];
        if i==1
            rotationMatrix1=rotationVectorToMatrix(vector1);
            axis2(i,j,1:3)=a;
        else
            axis2(i,j,1:3)=rotationMatrix1*a';  %��ǰ����a������һ������ax_rot��ת-phR(i,j)/2�Ƕȡ�
            Qa=axis2(i,j,1);    Ua=axis2(i,j,2);   Va=axis2(i,j,3);
            d1= -phR(i,j)/2.0*[Qa,Ua,Va]; %[Qa,Ua,Va]����ת��Ĺ��ᣬ
            rotationMatrix1 = rotationVectorToMatrix(d1)*rotationMatrix1;
%                        rotationMatrix1 = rotationVectorToMatrix(d1);
            % post-process to visualization purpose
            Q2=Qa/sqrt(Qa^2+Ua^2+Va^2);  U2=Ua/sqrt(Qa^2+Ua^2+Va^2);  V2=Va/sqrt(Qa^2+Ua^2+Va^2); 
            V2(isnan(V2)) = 0;
%              Q2=Qa/sqrt(Qa^2+Ua^2);  U2=Ua/sqrt(Qa^2+Ua^2);   
%              Q2(isnan(Q2)) = 0;U2(isnan(U2)) = 0;
            axiss2(i,j,:)=cat(3,Q2,U2,V2);
        end
    end
end

% bfloaxis3D=(axiss2+1)/2;
bfloaxis3D=axiss2;
bfphrr=phR;
bfloaxis3D(isnan(bfloaxis3D)) = 0;   
bfphrr(isnan(bfphrr)) = 0;   
% if showimg
%       figure, imshow(Stru,[]);
%       figure;    imagesc(phR);
%       figure;    imagesc((loaxis2+1)/2); %loaxis2 
%     figure;    imagesc((loaxis3+1)/2); %loaxis2 
%     figure;    imagesc((axiss2+1)/2); % axiss2
% end
% PhRcmap = jet(256);
% PRR = uint8(ind2rgb(uint8(phR/pi/0.5*256),PhRcmap)*256);
% bfPR=PRR;
LA=bfloaxis3D*0;PhR=bfphrr*0;
for j=1:size(Stru,2)
    LA(test_seg_top(j):size(LA,1),j,:)=bfloaxis3D(1:size(LA,1)-test_seg_top(j)+1,j,:);
    PhR(test_seg_top(j):size(PhR,1),j)=bfphrr(1:size(PhR,1)-test_seg_top(j)+1,j);
end
end