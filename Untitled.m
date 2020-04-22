clc;
clear all;
close all;
tic;
%% ����
pixel=0.2;  %%4.68um����    ��Ӧ͸������Ϊ4.5mm  961.5pixel
ratio=2/pixel;
range = 2.001*ceil(3273*ratio)+1;  %%708�����Ϊ���⻷�ߴ� 2.05�ɱ�
row=round(range)+1;
column=round(range)+1;
offset=(row+1)/2;
w=100;
lamda=488e-3;
k=2*pi/lamda;
%% ����mask
ring=zeros(row,column);%Ԥ����ƽ��I64�ĻҶ�ֵΪ0
[m,n]=meshgrid(linspace(-row/2,row/2-1,row));%ȷ������ϵ������ԭ���λ��
r3=ceil(2802*ratio);
% r3=ceil(490*ratio);
% r4=ceil(500*ratio);
r4=ceil(3273*ratio);
a=0;b=0;%Բ��λ�ÿ��Ʊ���a��Բ��λ�ÿ��Ʊ���b
D=((m+a).^2+(n+b).^2).^(1/2);%Բ������ϵʽ
i=find((D>=r3)&(D<=r4));%f������������D<=r1,D>=r2�����ص�ĵ��±�ֵ
ring(i)=1;

% 
% r1=ceil(198*ratio);
% r2=ceil(166*ratio);%Բ�İ뾶��С���Ʊ���r(��λPixel)
% a=0;b=0;%Բ��λ�ÿ��Ʊ���a��Բ��λ�ÿ��Ʊ���b
% D=((m+a).^2+(n+b).^2).^(1/2);%Բ������ϵʽ
% i=find((D<=r1)&(D>=r2));%f������������D<=r1,D>=r2�����ص�ĵ��±�ֵ
% ring(i)=1;%���ص㸳ֵ (1:N,1:N)
% 
% rect((column/2-row/2+1):(column/2+row/2),(column/2-row/2+1):(column/2+row/2))=ring;
% ring(:,1:(offset-w/2-1))= 0;
% ring(:,(offset+w/2+1):row)= 0;
% 
% ring(12230:12240,:)= 0;%ȥ������Բ��
% ring(1:11,:)= 0;
figure
imshow(ring);
title('figure1:mask');
%% 

% ��������
L0=(row-1)*pixel;                           %��ֵ������ߴ磬��λ:΢��
x=linspace(-L0/2,L0/2,row);                 %���ɱ��ξ� x������
y=linspace(-L0/2,L0/2,row);                 %���ɱ��ξ� y������
[x,y]=meshgrid(x,y);                        %���ξ�����ʵ�ߴ磬��λum

% ����
% r0=6.108e3; %waist radius selected��Chameloen�Ĺ�߳�ʼֱ����1.2mm�������󸲸�SLM���ȫ��1080x1080���أ��Խ���ֱ��12.217mm��6.108/0.008=763.5
% A0=(1./(pi*r0^2)).*exp(-((x).^2+(y).^2)./(r0^2));%Gauss beam
A0=1;                                   %%ƽ�й���
[xx,yy]=meshgrid(-(row-1)/2:(row-1)/2,-(row-1)/2:(row-1)/2);
Aperture=(1-sign(xx.^2+yy.^2-((row-1)./2).^2))./2.*abs(sign(xx.^2+yy.^2-((row-1)./2).^2));
A0=A0.*Aperture;
figure;
imshow(A0);
title('figure2:incident beam');

% ������Ĥ�����䵽͸��ǰ����
g=A0.*ring;
d=15e3;                                %��ֵ�۲�����������ľ���,��λ:΢�� 
kethi=linspace(-1./2./L0,1./2./L0,row).*row; %����Ƶ������
nenta=linspace(-1./2./L0,1./2./L0,column).*column;
[kethi,nenta]=meshgrid(kethi,nenta);
H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %���������䴫�ݺ���H
fa=fftshift(fft2(g));                 %�������Ϲⳡ�ĸ���Ҷ�任
Fuf=fa.*H;                            %�ⳡ��Ƶ���봫�ݺ������
U=ifft2(Fuf);                         %�ڹ۲����ϵĹⳡ�ֲ�
Am=abs(U).^2;                            %�ڹ۲����ϵĹ�ǿ�ֲ�
figure;
surf(x,y,Am), axis equal, axis tight, view([0, 270]), colorbar,colormap('hot'), shading interp;
title('figure4:in front of fourier lens');

% pattern_before_lens_line=Am(3000,round(offset)-3001:round(offset)+2990);
% figure
% plot(pattern_b efore_lens_line);
% title('figure8��pattern');

% U(:,1:(offset-250-1))= 0;
% U(:,(offset+250+1):row)= 0;
% imshow(U)
% % ����͸�������䵽����
f=15e3; %͸������100mm
g=U.*exp(-1i.*k./2./f.*(x.^2+y.^2)); %����͸����͸���ʺ���

% ��Z�����䴫��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��Z�����䴫��
 kethi=linspace(-1./2./L0,1./2./L0,row).*row; %����Ƶ������
 nenta=linspace(-1./2./L0,1./2./L0,column).*column;
 [kethi,nenta]=meshgrid(kethi,nenta);
 ly=20;%z����������
 zstart=14.98e3;
 zfinish=15.02e3;%%%%z�������ʼλ�ú���ֹλ��
 h=waitbar(0,'Cross section Caculating ...');
 index=1;
  fa=fft2(g);              %�������Ϲⳡ�ĸ���Ҷ�任
for d=linspace(zstart,zfinish,ly)
    H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %���������䴫�ݺ���H

    Fuf=fa.*H;                         %�ⳡ��Ƶ���봫�ݺ������
    U=ifft2(Fuf);                      %�ڹ۲����ϵĹⳡ�ֲ�
    Am=abs(U).^2;                         %�ڹ۲����ϵĹ�ǿ�ֲ�
    A_z(:,index)=diag(Am,0);
    A_z(:,index)=Am(ceil(row/2),:);
    A_y(:,index)=diag(Am,0);
    A_y(:,index)=Am(:,ceil(column/2));
    index=index+1;
    waitbar((index-1)./ly);
end
close(h);
[z2,x2]=meshgrid(zstart:((zfinish-zstart)./(ly-1)):zfinish,linspace(-L0/2,L0/2,row));
figure
surf(z2,x2,A_z),  colormap('hot'),shading interp; title('figure5:crelay range');

[z2,y2]=meshgrid(zstart:((zfinish-zstart)./(ly-1)):zfinish,linspace(-L0/2,L0/2,row));
figure
surf(z2,y2,A_y),  colormap('hot'),shading interp; title('figure5:zrelay range');



x_z=A_z(round(offset),:);
figure
plot(x_z);
title('figure6��c����������ֲ�������Χ');
% 
x_y=A_y(round(offset),:);
figure
plot(x_y);
title('figure6��z����������ֲ�������Χ');

bessel_line=A_y(round(offset)-101:round(offset)+99,10);
figure
plot(bessel_line);
title('figure8�����������ȫ��');

%% ����������ֲ�
d=15e3;
kethi=linspace(-1./2./L0,1./2./L0,row).*row; %����Ƶ������
nenta=linspace(-1./2./L0,1./2./L0,column).*column;
[kethi,nenta]=meshgrid(kethi,nenta);
H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %���������䴫�ݺ���H
fa=fft2(g);                           %�������Ϲⳡ�ĸ���Ҷ�任
Fuf=fa.*H;                            %�ⳡ��Ƶ���봫�ݺ������
U=ifft2(Fuf);                         %�ڹ۲����ϵĹⳡ�ֲ�
Am=abs(U).^2;                        %�ڹ۲����ϵĹ�ǿ�ֲ�



o = linspace(-500,500,1001);    
p = linspace(-500,500,1001); 
figure;
A = Am(round(offset)-501:round(offset)+499,round(offset)-501:round(offset)+499);
surf(o,p,A), 
axis equal;
axis tight, 
view([0, 270]),
colorbar,colormap('hot'), shading interp;title('figure7:bessel');

bessel_line=Am(round(offset)-101:round(offset)+99,round(offset));
figure
plot(bessel_line);
title('figure8�����������ȫ��');

bessel_line=Am(round(offset),round(offset)-101:round(offset)+99);
figure
plot(bessel_line);
title('figure8�������������');


% %% �����ڶ�ö͸��
% f=125e3; %͸������100mm
% g=U.*exp(-1i.*k./2./f.*(x.^2+y.^2)); %����͸����͸���ʺ���
% d=125e3;
% kethi=linspace(-1./2./L0,1./2./L0,row).*row; %����Ƶ������
% nenta=linspace(-1./2./L0,1./2./L0,column).*column;
% [kethi,nenta]=meshgrid(kethi,nenta);
% H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %���������䴫�ݺ���H
% fa=fft2(g);                           %�������Ϲⳡ�ĸ���Ҷ�任
% Fuf=fa.*H;                            %�ⳡ��Ƶ���봫�ݺ������
% U=ifft2(Fuf); 
% 
% %% �������������ﾵ
% f=10e3; %͸������100mm
% g=U.*exp(-1i.*k./2./f.*(x.^2+y.^2)); %����͸����͸���ʺ���
% d=10e3;
% kethi=linspace(-1./2./L0,1./2./L0,row).*row; %����Ƶ������
% nenta=linspace(-1./2./L0,1./2./L0,column).*column;
% [kethi,nenta]=meshgrid(kethi,nenta);
% H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %���������䴫�ݺ���H
% fa=fft2(g);                           %�������Ϲⳡ�ĸ���Ҷ�任
% Fuf=fa.*H;                            %�ⳡ��Ƶ���봫�ݺ������
% U=ifft2(Fuf);                         %�ڹ۲����ϵĹⳡ�ֲ�
% Am=abs(U).^2;  
% figure
% plot(Am(:,round(offset)));
% title('figure6:bessel_plot');


figure;
surf(x,y,Am), 
axis equal;
axis tight, 
view([0, 270]),
colorbar,colormap('hot'), shading interp;title('bessel');