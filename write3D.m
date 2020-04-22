function write3D(image, name)
%WRITE3D writes given 3D image to given file name
%   Inputs:
%       image - 3D image stack
%       name - file name string
    imwrite(uint16(image(:,:,1)), name, 'tiff');
    for i = 2:size(image, 3)
        imwrite(uint16(image(:,:,i)), name, 'tiff', 'WriteMode', 'append');
    end
end

