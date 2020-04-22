clc;
clear all;
close all;
for i1=1:1
    for i2=4:4
        for i3=2:2
            a='raw';
            b='.tif';
            filename2 = strcat(a,num2str(i1),num2str(i2),num2str(i3),b);
            testpara(filename2);
        end
    end
end
    