function mainCS3D(filename,beta, meanStd, xySize, zSize)
clc
 disp(['Now is calculating ',filename,' and beta is set as ',num2str(beta)]);
%% IMAGE INPUT
% Choose parameters;
    XRes = xySize;
    YRes = xySize;
    ZRes = zSize;

if meanStd > 0.08
    gamma = 0.5;
    blurSigma = 1e-3; % no Gaussian blur
    xySensor = 2;
    zSensor = 2;
else
    gamma = 1;
    blurSigma = 1;
    xySensor = 3.25;
    zSensor = 4;
end
inputfile = zeros(XRes, YRes, ZRes);
for i = 1:ZRes
    im1=imread(filename, i);
    %inputimg=rgb2gray(im1); 
    im1=im2double(im1);
    inputfile(:,:,i) =im1;
end

inputfile = inputfile./max(max(max(inputfile)));

%% PROCESSING PARAMETERS
FactorX =4;
FactorY = FactorX;
FactorZ =4;



HRX = FactorX * XRes;
HRY = FactorY * YRes; 
HRZ = FactorZ * ZRes;

%% CREATE INTERPOLATED GRID
img_interp = resize3D(inputfile, HRX, HRY, HRZ);

% %% GENERATE PSF
% psf = BesFilt3D(xySensor, zSensor/2, 1.1,2, FactorX, FactorZ);
psf = BesFilt3D(xySensor, zSensor/2, xySensor*2, zSensor*1.2, FactorX, FactorZ);
% psf_name ='psf_imaged_2.tif';
% psf = zeros(15,15,50);
% for i = 1:50
%     psf_input=imread(psf_name, i);
%     psf_input=im2double(psf_input);
%     psf(:,:,i) =psf_input;
% end

%psf=Bessepsf();
psf = modpsf(psf);
psf = psf ./ sum(sum(sum(psf)));



     
n= [HRX, HRY, HRZ];
A  = csConv(n,psf); % A
At = A';                % transpose of A

%% RUN OPTIMIZATION
target = img_interp;
lambda = find_lambdamax_l1_ls_nonneg(At,target(:));



[img_est, status] = l1_ls_nonneg(A,At,HRX*HRY*HRZ, HRX*HRY*HRZ, target(:), beta*lambda, 0.01);
% [img_est, status] = l1_ls_nonneg(A,At,HRX*HRY*HRZ, HRX*HRY*HRZ, target(:), 0.8, 0.01);
gather(img_est);
img_est = reshape(img_est, [(size(img_interp,1)) (size(img_interp,2))  (size(img_interp,3))]);



img_est = img_est .^ gamma;
img_est = imgaussfilt3(img_est, blurSigma);

%figure, imshow(img_est,[0,1]);

output16 = uint16(img_est.*(65535/max(img_est(:))/4));
output8 = uint8(img_est.*((255/max(img_est(:)))/4));


write3D(output16,[filename,'_X',num2str(FactorX),'Z',num2str(FactorZ),'interp_',num2str(beta),'lambda16bit0.001_Bespsf.tif']);
%write3D(output8,'raw_4xinterp_0.5lambda8bit.tif');
%write3D(img_est,'raw_4xinterp_0.5lambdadouble.tif');



   

