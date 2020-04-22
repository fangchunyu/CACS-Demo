

filename = '3D.tif';
XRes = 10;
YRes = 10;
ZRes = 15;

%90 90 80

%% PROCESSING PARAMETERS
FactorX = 3;
FactorY = FactorX;
FactorZ = 2;

xySensor = 1.6;
zSensor = 6;

HRX = FactorX * XRes;
HRY = FactorY * YRes;
HRZ = FactorZ * ZRes;


%% GENERATE PSF
psf = gausFilt3D(xySensor, zSensor/2, xySensor*2, zSensor*1.2, FactorX, FactorZ);
psf = modpsf(psf);
psf = psf ./ sum(sum(sum(psf)));
start1 = (size(psf,1)+1)/2;
start2 = (size(psf,2)+1)/2;
start3 = (size(psf,3)+1)/2;


%%
ball1 = [10 9 10 2];
ball2 = [10 15 10 2];
sq1 = [10 20 20 1 5 1];

rawimage = zeros(HRX, HRY, HRZ);
for i=1:HRX
    for j=1:HRY
        for k=1:HRZ
            %ball1
            dist = sqrt((i-ball1(1))^2 + (j-ball1(2))^2 + (k-ball1(3))^2);
            if(  dist <= ball1(4))
                rawimage(i,j,k) =200-dist*10;
            end
            
            %ball2
            dist = sqrt((i-ball2(1))^2 + (j-ball2(2))^2 + (k-ball2(3))^2);
            if(  dist <= ball2(4))
                rawimage(i,j,k) =200-dist*10;
            end
            
            %sq1
            if( (abs(i-sq1(1))<=sq1(4)) && (abs(j-sq1(2))<=sq1(5))  && (abs(k-sq1(3))<=sq1(6)) )
                rawimage(i,j,k) =160;
            end
        end
    end
end


write3D(rawimage,'rawimage.tif');

%% Convolve
convimage = convn(rawimage, psf);
hrimage = convimage( start1:start1+HRX-1,start2:start2+HRY-1,start3:start3+HRZ-1);

write3D(hrimage,'blurrimage.tif');
%% downsample
lrimage = zeros(XRes, YRes, ZRes);
for i=1:XRes
    for j=1:YRes
        for k=1:ZRes
            lrimage(i,j,k) = sum(sum(sum(hrimage((i-1)*FactorX+1:i*FactorX,(j-1)*FactorY+1:j*FactorY,(k-1)*FactorZ+1:k*FactorZ))));
            lrimage(i,j,k) =    lrimage(i,j,k)/(FactorX*FactorY*FactorZ);
        end
    end
end
lrimage = lrimage./(max(max(max(lrimage))));
lrimage = lrimage.*65535;
write3D(lrimage,filename);
img_interp = resize3D(lrimage, HRX, HRY, HRZ);