clc;
clear all;
close all;



% 
 for i=1:1
%     
 filename = strcat('D:\example_data\line_like_neurons\LR\thy1dense.tif');
 
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
% imgBlock =65535*ones(50,50,50);
% imgB=uint8(imgBlock./256)
% val=entropy(imgB);
img = sort(imgBlock(:), 'descend');
signal = img(1: 0.01*xySize*xySize*zSize);
StdDev = std(signal);
meanStd = mean(StdDev(:));
mainCS3D(filename, beta, meanStd, xySize, zSize);
    % TODO： 保存推测最佳lambda下的图像结果 img
    
%     for j=1:10
%         lambda = j/10;
%         img = CS3D(imgBlock, lambda);
%         % TODO： 保存不同lambda下的图像结果 img
% %     end        
 end