function  res = csConv(n,psf)

res.adjoint = 0;
res.n = n;
res.psf = psf;

rpsf = zeros(size(psf));
if(length(n)==2)
    for i=1:size(psf,1)
        for j= 1: size(psf,2)
            rpsf(i,j) = psf(size(psf,1)+1-i, size(psf,2)+1-j);
        end
    end
end

if(length(n)==3)
    for i=1:size(psf,1)
        for j= 1: size(psf,2)
            for k=1:size(psf,3)
                rpsf(i,j, k) = psf(size(psf,1)+1-i, size(psf,2)+1-j, size(psf,3)+1-k);
            end
        end
    end
end
res.rpsf = rpsf;

res = class(res,'csConv');
