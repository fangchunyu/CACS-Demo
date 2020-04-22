function newpsf = modpsf(psf)
dim = length(size(psf));

psfold = double(psf);

if(dim ==2)
    if( (mod(size(psfold,1),2)==0)&&(mod(size(psfold,2),2)==0))
        psf = zeros(size(psfold,1)-1, size(psfold,2)-1);
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                psf(i,j) = (psfold(i,j) + psfold(i+1,j) + psfold(i,j+1) + psfold(i+1,j+1))/4;
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==0)&&(mod(size(psfold,2),2)==1))
        psf = zeros(size(psfold,1)-1, size(psfold,2));
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                psf(i,j) = (psfold(i,j) + psfold(i+1,j))/2;
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==1)&&(mod(size(psfold,2),2)==0))
        psf = zeros(size(psfold,1), size(psfold,2)-1);
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                psf(i,j) = (psfold(i,j) + psfold(i,j+1))/2;
            end
        end
    end
end



if(dim ==3)
    if( (mod(size(psfold,1),2)==0)&&(mod(size(psfold,2),2)==0)&&(mod(size(psfold,3),2)==0))
        psf = zeros(size(psfold,1)-1, size(psfold,2)-1, size(psfold,3)-1);
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j) = (psfold(i,j,k) + psfold(i+1,j, k) + psfold(i,j+1, k)  + psfold(i,j,k+1) + psfold(i+1,j+1,k) + psfold(i+1,j,k+1) + psfold(i,j+1,k+1)+ psfold(i+1,j+1, k+1))/8;
                end
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==0)&&(mod(size(psfold,2),2)==1)&&(mod(size(psfold,3),2)==1))
        psf = zeros(size(psfold,1)-1, size(psfold,2),  size(psfold,3));
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j, k) = (psfold(i,j, k) + psfold(i+1,j, k))/2;
                end
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==1)&&(mod(size(psfold,2),2)==0)&&(mod(size(psfold,3),2)==1))
        psf = zeros(size(psfold,1), size(psfold,2)-1,  size(psfold,3));
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j, k) = (psfold(i,j,k) + psfold(i,j+1, k))/2;
                end
            end
        end
    end
    if( (mod(size(psfold,1),2)==1)&&(mod(size(psfold,2),2)==1)&&(mod(size(psfold,3),2)==0))
        psf = zeros(size(psfold,1), size(psfold,2),  size(psfold,3)-1);
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j, k) = (psfold(i,j,k) + psfold(i,j, k+1))/2;
                end
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==0)&&(mod(size(psfold,2),2)==0)&&(mod(size(psfold,3),2)==1))
        psf = zeros(size(psfold,1)-1, size(psfold,2)-1,  size(psfold,3));
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j, k) = (psfold(i,j, k) + psfold(i+1,j, k)+ psfold(i,j+1, k) + psfold(i+1,j+1, k) )/4;
                end
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==0)&&(mod(size(psfold,2),2)==1)&&(mod(size(psfold,3),2)==0))
        psf = zeros(size(psfold,1)-1, size(psfold,2),  size(psfold,3)-1);
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j, k) = (psfold(i,j, k) + psfold(i+1,j, k)+ psfold(i,j, k+1) + psfold(i+1,j, k+1) )/4;
                end
            end
        end
    end
    
    if( (mod(size(psfold,1),2)==1)&&(mod(size(psfold,2),2)==0)&&(mod(size(psfold,3),2)==0))
        psf = zeros(size(psfold,1), size(psfold,2)-1,  size(psfold,3)-1);
        for i=1:size(psf,1)
            for j= 1:size(psf,2)
                for k=1:size(psf,3)
                    psf(i,j, k) = (psfold(i,j, k) + psfold(i,j+1, k)+ psfold(i,j, k+1) + psfold(i,j+1, k+1) )/4;
                end
            end
        end
    end
    
end


newpsf = psf;