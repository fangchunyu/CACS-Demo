function psf=Bessepsf()
filename='bessel psf.tif';
XRes = 7;
YRes =7;
ZRes = 7;
inputfile = zeros(XRes, YRes, ZRes);
for i = 1:ZRes
  
    im1=imread(filename, i);
    %inputimg=rgb2gray(im1); 
    im1=im2double(im1);
    inputfile(:,:,i) =im1;
end

X_H=14;
Y_H=14;
Z_H=17;
[x,y,z]=meshgrid(1:XRes,1:YRes,1:ZRes);
d_xy=XRes/X_H;
d_z=ZRes/Z_H;
[x1,y1,z1]=meshgrid(1:d_xy:XRes,1:d_xy:YRes,1:d_z:ZRes);
psf=interp3(x,y,z,inputfile,x1,y1,z1);
% psf=zeros(13,13,15);
% psf(4:10,4:10,5:11)=inputfile(1:7,1:7,1:7);
