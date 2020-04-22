function lambda = CalcLambda32(imgblock)
% For line_like
% nhoodsz = 5;

% For point_like
nhoodsz = 3;
nhood = ones(nhoodsz,nhoodsz,nhoodsz);
test1 = entropyfilt(uint8(imgblock./256),nhood);


aa = reshape(test1,[],1);
c = sort(aa,'descend');
cc = c(1:5);
val = mean(cc);

% For line_like
lambda = (-0.0841)*val + 0.5199;

% For point_like
% lambda = (-0.3538)*val + 1.8224


end