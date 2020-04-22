function lambda = CalcLambda32(imgblock)

% nhoodsz = 15;
nhoodsz = 3;
nhood = ones(nhoodsz,nhoodsz,nhoodsz);
test1 = entropyfilt(uint8(imgblock./256),nhood);


aa = reshape(test1,[],1);
c = sort(aa,'descend');
cc = c(1:5);
val = mean(cc);
% val=c;
% lambda = (-0.085)*val + 0.5199;
lambda = (-0.3538)*val + 1.8224;
% lambda = (-0.14)*val + 0.71;

end