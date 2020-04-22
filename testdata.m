clc;
clear all;
close all;
FactorX =2;
FactorY = FactorX;
FactorZ =2;

xySensor = 1.032;
zSensor = 1;

% HRX = FactorX * XRes;
% HRY = FactorY * YRes; 
% HRZ = FactorZ * ZRes;

%% CREATE INTERPOLATED GRID
%img_interp = resize3D(inputfile, HRX, HRY, HRZ);

%% GENERATE PSF
psf = BesFilt3D(xySensor, zSensor/2, xySensor*2, zSensor*1.2, FactorX, FactorZ);