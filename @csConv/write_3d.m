function write_3d(result,name,Z)

maxres = max(max(max(result)));
result = result./maxres;
result = result.*65535;

%save tiff
for zz = 1:Z
    Imgout = uint16(result(:,:,zz));
    if zz==1
        imwrite(Imgout, name);
    else
        imwrite(Imgout, name, 'WriteMode','append');
    end
end
