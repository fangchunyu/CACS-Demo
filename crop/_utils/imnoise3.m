function I = imnoise3(stack,variance)

for k = 1 : size(stack, 3)
    stack(:,:,k) = imnoise(stack(:,:,k),'gaussian', 0, variance);
end

I = stack;