clc;
XRes =210;
YRes =210;
ZRes =210;
% inputfile = zeros(XRes, YRes, ZRes);
% for i = 1:ZRes
%   
%     im1=imread(filename, i);
%     %inputimg=rgb2gray(im1); 
%     im1=im2double(im1);
%     inputfile(:,:,i) =im1;
% end
% 
% inputfile = inputfile./max(max(max(inputfile)));

%% PROCESSING PARAMETERS
FactorX = 2;
FactorY = FactorX;
FactorZ =2;

xySensor = 1.032;
zSensor = 1;

HRX = FactorX * XRes;
HRY = FactorY * YRes; 
HRZ = FactorZ * ZRes;

%% CREATE INTERPOLATED GRID
%img_interp = resize3D(inputfile, HRX, HRY, HRZ);

%% GENERATE PSF
psf = gausFilt3D(xySensor, zSensor/2, xySensor*2, zSensor*1.2, FactorX, FactorZ);
psf = modpsf(psf);
psf = psf ./ sum(sum(sum(psf)));


     
% n= [HRX, HRY, HRZ];
% A  = csConv(n,psf); % A
% At = A';                % transpose of A