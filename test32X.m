% Test 3.2X

imgSize = 100;
factor = 4;
blockSize = 25;
testTimes = 10; %����10��

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
    % TODO�� ����ԭʼͼ��� imgLargeBlock
    
    imgBlock = imgraw(rx:(rx+blockSize-1),ry:(ry+blockSize-1),rz:(rz+blockSize-1));
    %imgLargeBlock��imgBlock�ֱ���12.6X��3.2X�������Ӧ��������
    
    %�������lambda
    lambda = CalcLambda32(imgBlock);
%     img = CS3D(imgBlock, lambda);

 write3D(imgBlock,['RAW',num2str(i),' 3.2X.tif']);
    write3D(imgLargeBlock,['RAW',num2str(i),' 12.6X.tif']);
    filename = strcat('RAW',num2str(i),' 3.2X.tif');
    mainCS3D(filename, lambda);
    % TODO�� �����Ʋ����lambda�µ�ͼ���� img
    
    for j=1:9
        lambda = j/10;
        mainCS3D(filename, lambda);
        % TODO�� ���治ͬlambda�µ�ͼ���� img
    end        
end