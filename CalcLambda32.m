function lambda = CalcLambda32(nh, k1, b1, imgblock)


nhoodsz = nh;
nhood = ones(nhoodsz,nhoodsz,nhoodsz);
% nhoodsz = 15;


test1 = entropyfilt(uint8(imgblock./256),nhood);


aa = reshape(test1,[],1);
c = sort(aa,'descend');
cc = c(1:5);
val = mean(cc);
% val=c;
lambda = k1 * val + b1;


end