function h = BesFilt3D(xySensor, zSensor, xyReal, zReal, xyImprove, zImprove)
%
%gausFilt3D Generates a 3D gaussian filter
%   Inputs:
%   xySensor - X and Y pixel resolution
%   zSensor - Z pixel resolution
%   xyReal - X and Y optical resolution
%   zReal - Z optical resolution
%   xyImprove - X and Y improvement factor
%   zImprove - Z improvement factor
%   Outputs:
%   h - 3D gaussian filter


%% Compute filter size
Size = 100; % Maximum X,Y,Z pixel size
xyGrid = round(xyReal/(xySensor/xyImprove));
zGrid = round(zReal/(zSensor/zImprove));
xySigma = round(xyGrid/2)-1;
% xySigma =6;

if xySigma < 1
    xySigma = 1;
end
zSigma = round(zGrid/2)-1;
% zSigma = 6;
if zSigma < 1
    zSigma = 1;
end

%% Generate filter slices
xyH = fspecial('gaussian', [Size Size], xySigma);
xyH = xyH * (65355/max(xyH(:)));
zHLine = fspecial('gaussian', [1 Size], zSigma);



% zline = besselj(0, 1 : floor(Size/2 ));
% 
% zHLine = [flip(zline) zline];
% zHLine(zHLine < 0) = 0;
% 
 zHLine = zHLine * (65535/max(zHLine(:)));

h = zeros(Size, Size, Size);
for i = 1:Size
    h(:,:,i) = xyH * zHLine(1,i);
end

h = h * (65355/max(h(:)));

%% Crop filter to 1.5x size
xyStart = round(Size/2-xyGrid*1.5);
xyEnd = xyStart+(xyGrid*3);
zStart = round(Size/2-zGrid*1.5);
zEnd = zStart+(zGrid*3);
h = h(xyStart:xyEnd, xyStart:xyEnd, zStart:zEnd);

imwrite(uint16(h(:,:,1)'), 'psf1.tif', 'tiff');
for i = 2:size(h, 3)
    imwrite(uint16(h(:,:,i)'), 'psf1.tif', 'tiff', 'WriteMode', 'append');
end

end