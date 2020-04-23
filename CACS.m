clc;
clear all;
close all;



% 
 for i=1:1
%     
 filename = strcat('D:\example_data\point_like_cell_nuclei\LR\pidense.tif');
 
 if contains(filename, 'thy1')
     blockSize = 50;
     nh = 5;
     k1 = -0.085;
     b1 = 0.5199;
 else
     blockSize = 10;
     nh = 3;
     k1 = -0.3538;
     b1 = 1.8224;
 end
%  filename = 'test2.tif';
imgBlock = imgread3D(filename, blockSize); 
beta = CalcLambda32(nh, k1, b1, imgBlock);
imgBlock = imgBlock./max(max(max(imgBlock)));
xySize = size(imgBlock, 1);
zSize = size(imgBlock, 3);
img = sort(imgBlock(:), 'descend');
signal = img(1: 0.01*xySize*xySize*zSize);
StdDev = std(signal);
meanStd = mean(StdDev(:));
mainCS3D(filename, beta, meanStd, xySize, zSize);
      
 end