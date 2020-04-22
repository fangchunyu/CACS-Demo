%function testpara(filename1)
clc;
clear all;
close all;
tic
% filename2='bessel-2_deblurr.tif';
beta=[0.05,0.1,0.15];
%   for n=0.5
%       filename1 =  strcat('D:\CS_PRO\SI6\r',num2str(n),'.tif');
  filename1='\\192.168.1.11\d\besselcs\20190820\2\roi\6\488_lr_6.tif';
    for i=1:3
    mainCS3D(filename1,beta(i));
    end
%  end
% % % for i=1:5
% %   mainCS3D(filename1,beta(i));
%  end
toc

% filename1='Center.tif';
% 
% beta=[0.52,0.54,0.56,0.58,0.62,0.64,0.66,0.68,0.7,0.8];
% 
% % filename2='bessel-2_deblurr.tif'
% for i=1:10
%   mainCS3D(filename1,beta(i));
% end
% 
% filename1='Right.tif';
% % filename2='bessel-2_deblurr.tif';
% beta=[0.36,0.38,0.42,0.44,0.46,0.48,0.52,0.54,0.7,0.8];
% 
% for i=1:10
%   mainCS3D(filename1,beta(i));
% end

