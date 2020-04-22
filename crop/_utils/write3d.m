function write3d(img, name, bitdepth)
%WRITE3D writes given 3D image to given file name
%   Inputs:
%       image - 3D image stack
%       name - file name string
%       bitdepth - 8 or 16

%   Temporary : image2wt - the matrix used to write a tif image,in case
%   that the original image is changed.


%assert(max(img(:)) <= 255, 'image data to be saved must be in 8-bit')
if (bitdepth == 8)
    img = uint8(img);
elseif (bitdepth == 16)
    img = uint16(img);
else 
    error('unkonwn bitdepth: %d', bitdepth);
end
% image2wt = image;  
imwrite(img(:,:,1), name);

for i = 2:size(img, 3)
    imwrite((img(:,:,i)), name,  'WriteMode', 'append');
end

% for i = 1 : size(image2wt, 3)
%     imwrite(image2wt(:,:,i), sprintf('%d.tif', i));
% end


