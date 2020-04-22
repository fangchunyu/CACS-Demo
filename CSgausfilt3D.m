function h = CSgausfilt3D(xySensor, xyReal, xyImprove)
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
size = 100; % Maximum X,Y,Z pixel size
xyGrid = round(xyReal/(xySensor/xyImprove));
% zGrid = round(zReal/(zSensor/zImprove));
xySigma = round(xyGrid/2)-1;
if xySigma < 1
    xySigma = 1;
end
% zSigma = round(zGrid/2)-1;
% if zSigma < 1
%     zSigma = 1;
% end

%% Generate filter slices
xyH = fspecial('gaussian', [size size], xySigma);
xyH = xyH * (65355/max(xyH(:)));
% zHLine = fspecial('gaussian', [1 size], zSigma);
% zHLine = zHLine * (65355/max(zHLine(:)));

% h = zeros(size, size, size);
% for i = 1:size
%     h(:,:,i) = xyH * zHLine(1,i);
% end

xyH = xyH * (65355/max(xyH(:)));

%% Crop filter to 1.5x size
xyStart = round(size/2-xyGrid*1.5);
xyEnd = xyStart+(xyGrid*3);
% zStart = round(size/2-zGrid*1.5);
% zEnd = zStart+(zGrid*3);
xyH = xyH(xyStart:xyEnd, xyStart:xyEnd);
h = xyH;
end