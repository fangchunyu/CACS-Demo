%% resize XY of a 3-D stack
%  params:
%    -method: ['nearest', 'bicubic']

function I = resizeXY(stack, factor, method)


if length(size(stack)) ~= 3
    error('incorrect stack dimensions ');
end

slice = stack(:,:,1);
tmp = imresize(slice, factor, method);

nslices = size(stack, 3);
tmp = zeros([size(tmp), nslices]);

for k = 1 : nslices
    slice = stack(:,:,k);
    %tmp(:,:,k) = slice(1:k:end, 1:k:end);
    tmp(:,:,k) = imresize(slice, factor, method);
end

I = tmp;
    