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
% imgBlock =65535*ones(50,50,50);
% imgB=uint8(imgBlock./256)
% val=entropy(imgB);

beta = CalcLambda32(nh, k1, b1, imgBlock);
mainCS3D(filename, beta);
    % TODO： 保存推测最佳lambda下的图像结果 img
    
%     for j=1:10
%         lambda = j/10;
%         img = CS3D(imgBlock, lambda);
%         % TODO： 保存不同lambda下的图像结果 img
% %     end        
 end