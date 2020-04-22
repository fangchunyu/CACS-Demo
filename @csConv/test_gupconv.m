clc;
clear all;

brain_struct=load('brain25_100.mat');
brain=brain_struct.test;
psf=genpsf(4,1000);
result= GPUConv3D( single(brain), int32(size(brain)), single(psf), int32(size(psf)));
