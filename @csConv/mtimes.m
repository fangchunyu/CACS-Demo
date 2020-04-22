function res = mtimes(A,x)


%dimension of image
sz1 = A.n(1);
sz2 = A.n(2);
sz3 = A.n(3);
%dimension of psf
[psz1,psz2, psz3] = size(A.psf);
%note that size(psf,n) should be odd





if A.adjoint == 0 %A*x
    x = reshape(x,A.n);
    res = double(GPUConv3D( single(x), int32(size(x)), single(A.psf), int32(size(A.psf))));
else %At*x
    x = reshape(x,A.n);
    res = double(GPUConv3D( single(x), int32(size(x)), single(A.rpsf), int32(size(A.rpsf))));    
end


res = reshape(res, (sz1*sz2*sz3),[]);
