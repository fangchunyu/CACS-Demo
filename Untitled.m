clc;
clear all;
close all;
tic;
%% 参数
pixel=0.2;  %%4.68um像素    对应透镜焦距为4.5mm  961.5pixel
ratio=2/pixel;
range = 2.001*ceil(3273*ratio)+1;  %%708处替代为最外环尺寸 2.05可变
row=round(range)+1;
column=round(range)+1;
offset=(row+1)/2;
w=100;
lamda=488e-3;
k=2*pi/lamda;
%% 生成mask
ring=zeros(row,column);%预定义平面I64的灰度值为0
[m,n]=meshgrid(linspace(-row/2,row/2-1,row));%确定坐标系及坐标原点的位置
r3=ceil(2802*ratio);
% r3=ceil(490*ratio);
% r4=ceil(500*ratio);
r4=ceil(3273*ratio);
a=0;b=0;%圆心位置控制变量a，圆心位置控制变量b
D=((m+a).^2+(n+b).^2).^(1/2);%圆函数关系式
i=find((D>=r3)&(D<=r4));%f返回满足条件D<=r1,D>=r2的像素点的单下标值
ring(i)=1;

% 
% r1=ceil(198*ratio);
% r2=ceil(166*ratio);%圆的半径大小控制变量r(单位Pixel)
% a=0;b=0;%圆心位置控制变量a，圆心位置控制变量b
% D=((m+a).^2+(n+b).^2).^(1/2);%圆函数关系式
% i=find((D<=r1)&(D>=r2));%f返回满足条件D<=r1,D>=r2的像素点的单下标值
% ring(i)=1;%像素点赋值 (1:N,1:N)
% 
% rect((column/2-row/2+1):(column/2+row/2),(column/2-row/2+1):(column/2+row/2))=ring;
% ring(:,1:(offset-w/2-1))= 0;
% ring(:,(offset+w/2+1):row)= 0;
% 
% ring(12230:12240,:)= 0;%去掉顶部圆弧
% ring(1:11,:)= 0;
figure
imshow(ring);
title('figure1:mask');
%% 

% 生成坐标
L0=(row-1)*pixel;                           %赋值衍射面尺寸，单位:微米
x=linspace(-L0/2,L0/2,row);                 %生成变形镜 x轴坐标
y=linspace(-L0/2,L0/2,row);                 %生成变形镜 y轴坐标
[x,y]=meshgrid(x,y);                        %变形镜的真实尺寸，单位um

% 光束
% r0=6.108e3; %waist radius selected；Chameloen的光斑初始直径是1.2mm，扩束后覆盖SLM面板全部1080x1080像素，对角线直径12.217mm，6.108/0.008=763.5
% A0=(1./(pi*r0^2)).*exp(-((x).^2+(y).^2)./(r0^2));%Gauss beam
A0=1;                                   %%平行光束
[xx,yy]=meshgrid(-(row-1)/2:(row-1)/2,-(row-1)/2:(row-1)/2);
Aperture=(1-sign(xx.^2+yy.^2-((row-1)./2).^2))./2.*abs(sign(xx.^2+yy.^2-((row-1)./2).^2));
A0=A0.*Aperture;
figure;
imshow(A0);
title('figure2:incident beam');

% 经过掩膜，衍射到透镜前表面
g=A0.*ring;
d=15e3;                                %赋值观察屏到衍射面的距离,单位:微米 
kethi=linspace(-1./2./L0,1./2./L0,row).*row; %给出频域坐标
nenta=linspace(-1./2./L0,1./2./L0,column).*column;
[kethi,nenta]=meshgrid(kethi,nenta);
H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %菲涅尔衍射传递函数H
fa=fftshift(fft2(g));                 %衍射面上光场的傅里叶变换
Fuf=fa.*H;                            %光场的频谱与传递函数相乘
U=ifft2(Fuf);                         %在观察屏上的光场分布
Am=abs(U).^2;                            %在观察屏上的光强分布
figure;
surf(x,y,Am), axis equal, axis tight, view([0, 270]), colorbar,colormap('hot'), shading interp;
title('figure4:in front of fourier lens');

% pattern_before_lens_line=Am(3000,round(offset)-3001:round(offset)+2990);
% figure
% plot(pattern_b efore_lens_line);
% title('figure8：pattern');

% U(:,1:(offset-250-1))= 0;
% U(:,(offset+250+1):row)= 0;
% imshow(U)
% % 经过透镜，衍射到后焦面
f=15e3; %透镜焦距100mm
g=U.*exp(-1i.*k./2./f.*(x.^2+y.^2)); %乘以透镜的透过率函数

% 沿Z轴衍射传播
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%沿Z轴衍射传播
 kethi=linspace(-1./2./L0,1./2./L0,row).*row; %给出频域坐标
 nenta=linspace(-1./2./L0,1./2./L0,column).*column;
 [kethi,nenta]=meshgrid(kethi,nenta);
 ly=20;%z方向仿真层数
 zstart=14.98e3;
 zfinish=15.02e3;%%%%z仿真的起始位置和终止位置
 h=waitbar(0,'Cross section Caculating ...');
 index=1;
  fa=fft2(g);              %衍射面上光场的傅里叶变换
for d=linspace(zstart,zfinish,ly)
    H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %菲涅尔衍射传递函数H

    Fuf=fa.*H;                         %光场的频谱与传递函数相乘
    U=ifft2(Fuf);                      %在观察屏上的光场分布
    Am=abs(U).^2;                         %在观察屏上的光强分布
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
title('figure6：c贝塞尔轴向分布瑞丽范围');
% 
x_y=A_y(round(offset),:);
figure
plot(x_y);
title('figure6：z贝塞尔轴向分布瑞丽范围');

bessel_line=A_y(round(offset)-101:round(offset)+99,10);
figure
plot(bessel_line);
title('figure8：贝塞尔半高全宽');

%% 贝塞尔截面分布
d=15e3;
kethi=linspace(-1./2./L0,1./2./L0,row).*row; %给出频域坐标
nenta=linspace(-1./2./L0,1./2./L0,column).*column;
[kethi,nenta]=meshgrid(kethi,nenta);
H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %菲涅尔衍射传递函数H
fa=fft2(g);                           %衍射面上光场的傅里叶变换
Fuf=fa.*H;                            %光场的频谱与传递函数相乘
U=ifft2(Fuf);                         %在观察屏上的光场分布
Am=abs(U).^2;                        %在观察屏上的光强分布



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
title('figure8：贝塞尔半高全宽');

bessel_line=Am(round(offset),round(offset)-101:round(offset)+99);
figure
plot(bessel_line);
title('figure8：贝塞尔横截面');


% %% 经过第二枚透镜
% f=125e3; %透镜焦距100mm
% g=U.*exp(-1i.*k./2./f.*(x.^2+y.^2)); %乘以透镜的透过率函数
% d=125e3;
% kethi=linspace(-1./2./L0,1./2./L0,row).*row; %给出频域坐标
% nenta=linspace(-1./2./L0,1./2./L0,column).*column;
% [kethi,nenta]=meshgrid(kethi,nenta);
% H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %菲涅尔衍射传递函数H
% fa=fft2(g);                           %衍射面上光场的傅里叶变换
% Fuf=fa.*H;                            %光场的频谱与传递函数相乘
% U=ifft2(Fuf); 
% 
% %% 产生贝塞尔的物镜
% f=10e3; %透镜焦距100mm
% g=U.*exp(-1i.*k./2./f.*(x.^2+y.^2)); %乘以透镜的透过率函数
% d=10e3;
% kethi=linspace(-1./2./L0,1./2./L0,row).*row; %给出频域坐标
% nenta=linspace(-1./2./L0,1./2./L0,column).*column;
% [kethi,nenta]=meshgrid(kethi,nenta);
% H=exp(j*k*d.*(1-lamda.*lamda.*(kethi.*kethi+nenta.*nenta)./2)); %菲涅尔衍射传递函数H
% fa=fft2(g);                           %衍射面上光场的傅里叶变换
% Fuf=fa.*H;                            %光场的频谱与传递函数相乘
% U=ifft2(Fuf);                         %在观察屏上的光场分布
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