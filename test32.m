clc;
clear all;
close all;


blockSize = 10;
% 
 for i=1:1
%     
 filename = strcat('D:\paper_example\pi\25\sparseHIP_lr',num2str(i),'.tif');
%  filename = 'test2.tif';
imgBlock = imgread3D(filename, blockSize); 
% imgBlock =65535*ones(50,50,50);
% imgB=uint8(imgBlock./256)
% val=entropy(imgB);

lambda1 = CalcLambda32(imgBlock);
mainCS3D(filename, lambda1);
    % TODO�� �����Ʋ����lambda�µ�ͼ���� img
    
%     for j=1:10
%         lambda = j/10;
%         img = CS3D(imgBlock, lambda);
%         % TODO�� ���治ͬlambda�µ�ͼ���� img
% %     end        
 end