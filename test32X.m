% Test 3.2X

imgSize = 100;
factor = 4;
blockSize = 25;
testTimes = 10; %测试10次

imgraw = imgread3D('3.2X.tif', imgSize);
imgLarge = imgread3D('12.6X.tif',imgSize * factor);

for i=1:testTimes
    
    %Get a random point
    rx =  round((imgSize -blockSize-1)*rand()) +1;
    ry =  round((imgSize -blockSize-1)*rand()) +1;
    rz =  round((imgSize -blockSize-1)*rand()) +1;
    rxL = rx * factor;
    ryL = ry * factor;
    rzL = rz * factor;
    
    imgLargeBlock = imgLarge(rxL:(rxL+blockSize*factor-1),ryL:(ryL+blockSize*factor-1),rzL:(rzL+blockSize*factor-1));
    % TODO： 保存原始图像块 imgLargeBlock
    
    imgBlock = imgraw(rx:(rx+blockSize-1),ry:(ry+blockSize-1),rz:(rz+blockSize-1));
    %imgLargeBlock和imgBlock分别是12.6X和3.2X中随机对应的两个块
    
    %计算最佳lambda
    lambda = CalcLambda32(imgBlock);
%     img = CS3D(imgBlock, lambda);

 write3D(imgBlock,['RAW',num2str(i),' 3.2X.tif']);
    write3D(imgLargeBlock,['RAW',num2str(i),' 12.6X.tif']);
    filename = strcat('RAW',num2str(i),' 3.2X.tif');
    mainCS3D(filename, lambda);
    % TODO： 保存推测最佳lambda下的图像结果 img
    
    for j=1:9
        lambda = j/10;
        mainCS3D(filename, lambda);
        % TODO： 保存不同lambda下的图像结果 img
    end        
end