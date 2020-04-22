function lambda = CalcLambda64(imgblock)

% nhoodsz = 5;
nhoodsz = 5;
nhood = ones(nhoodsz,nhoodsz,nhoodsz);
test1 = entropyfilt(uint8(imgblock./256), nhood);


aa = reshape(test1,[],1);
c = sort(aa,'descend');
cc = c(1:50);
val = mean(cc);

lambda = (-0.17)*val + 0.82;

end