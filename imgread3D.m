function [img] = imgread3D(fname, z)

a = imread(fname);
[x y] = size(a);
b = zeros(x,y,z);
for i=1:z
    b(:,:,i) = imread(fname,i);
end
img = b;
end