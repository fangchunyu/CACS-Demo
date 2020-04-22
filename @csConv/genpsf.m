function psfout = genpsf(Sigma, cropThreshold)

gensize = 100; % Maximum X,Y,Z pixel size during generation
cropthershold = cropThreshold; % crop psf while pixel index is less than cropthershold (maximun = 65535)

xySigma = Sigma;
zSigma = xySigma;


%% Generate filter slices
xyH = fspecial('gaussian', [gensize gensize], xySigma);
xyH = xyH * (65535/max(xyH(:)));
zHLine = fspecial('gaussian', [1 gensize], zSigma);
zHLine = zHLine * (65535/max(zHLine(:)));

h = zeros(gensize, gensize, gensize);
for i = 1:gensize
    h(:,:,i) = xyH * zHLine(1,i);
end

h = h * (65535/max(h(:)));

%% Clip
psf = uint16(h);
sx = gensize;
sy = gensize;
sz = gensize;
psfm = psf(:,:,round(sz/2));
xcrop1 = 0;
for i=1:round(sx/2)
    if( max(psfm(i,:))<= cropthershold)
        xcrop1 = i;
    end
end
xcrop1 = xcrop1+1;

xcrop2 = gensize;
for i=round(sx/2):sx
    if( max(psfm(i,:)) > cropthershold)
        xcrop2 = i;
    end
end

xcrop1 = min(xcrop1, (gensize - xcrop2 +1));
xcrop2 = gensize - xcrop1 + 1;


ycrop1 = 0;
for i=1:round(sy/2)
    if( max(psfm(:,i))<= cropthershold)
        ycrop1 = i;
    end
end
ycrop1 = ycrop1+1;

ycrop2 = gensize;
for i=round(sy/2):sy
    if( max(psfm(:,i))> cropthershold)
        ycrop2 = i;
    end
end

ycrop1 = min(ycrop1, (gensize - ycrop2 +1));
ycrop2 = gensize - ycrop1 + 1;

output = psf(xcrop1:xcrop2, ycrop1:ycrop2, :);
output1 = output;

zcrop1 = 0;
for i=1:round(sz/2)
    if(max(max(output1(:,:,i)))<= cropthershold)
        zcrop1 = i;
    end
end
zcrop1 = zcrop1 +1;

zcrop = gensize;
for i=round(sz/2):sz
    if(max(max(output1(:,:,i)))> cropthershold)
        zcrop2 = i;
    end
end

zcrop1 = min(zcrop1, (gensize - zcrop2 + 1));
zcrop2 = gensize - zcrop1 + 1;

output = output1(:,:, zcrop1:zcrop2);
output = uint16(output);
outputold = output;


fp = fopen('./psf.raw','w');
ct1 = fwrite(fp,output, 'uint16');
fclose(fp);

temp = size(output);
image = output;

imwrite(uint16(image(:,:,1)'), 'psf.tif', 'tiff');
for i = 2:size(image, 3)
    imwrite(uint16(image(:,:,i)'), 'psf.tif', 'tiff', 'WriteMode', 'append');
end
imshow(output(:,:,temp(3)/2)',[]);

psfout = output;
