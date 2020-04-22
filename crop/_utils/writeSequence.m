function writeSequence(img, path)

assert(max(img(:)) <= 255, 'image data to be saved must be in 8-bit')
img = uint8(img);

if ~exist(path, 'dir')
    mkdir(path);
end
for i = 1 : size(img, 3)
    imwrite(img(:,:,i), [path, sprintf('\\%05d.tif', i)]);
end
end
