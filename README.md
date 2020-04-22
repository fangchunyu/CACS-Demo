# CACS-Demo
Content Aware Compressed Sensing for 3D resolution-enhancemeng of microscopy images 

## Requirements

CACS is built with Matlab. Technically there are no limits to the operation system to run the code, but Windows system is recommonded, on which the software has been tested.

The inference process of the CACS needs CUDA-enabled GPU device, which can speed up the inference. 

The inference process has been tested with:

 * Windows 10 pro (version 1903)
 * MatlabR2017a (64 bit)
 * Intel Xeon E5-2630 v4 CPU @2.80GHz
 * Nvidia GeForce GTX 1080 Ti

## Install

1. Install MatlabR2017a 
2. Install CUDA7.5.
3. Download the CACS-Demo.zip and unpack it. The directory tree should be: 

```  
CACS-Demo   
    .
    ├── @csConv
    ├── BesFilt3D.m
    ├── CalcLambda32.m
    ├── CSgausfilt3D.m
    ├── find_lambdamax_l1_ls_nonneg.m
    ├── gausFilt3D.m
    ├── example_data
        └── Line_like
        └── Point_like
```

## Usage

#### Inference

This toturial contains example data of line-like and point-like signals from Thy1-GFP-M and PI labelled mouse brain(see example_data/):

```
├── example_data
    └── Line_like
        └── LR
            └── thy1sparse.tif                (3.2x bessel-sheet ROI with sparse signals of Thy1-GFP-M mouse brain)
            └── thy1dense.tif                 (3.2x bessel-sheet ROI with dense signals of Thy1-GFP-M mouse brain)
        └── expected_outputs
    └── Point_like
        └── LR
            └── pisparse.tif                (3.2x bessel-sheet ROI with sparse signals of PI-labelled mouse brain)
            └── pidense.tif                 (3.2x bessel-sheet ROI with dense signals of PI-labelled mouse brain)
        └── expected_outputs

```

The expected outputs by the CACS of each input LR can be found in the corresponding 'expected_outputs' directory. 

To set parameters in the CACS, open these Matlab file :

```
CACS.m
CalcLambda32.m
mainCS3D.m

```

Setting parameters:

1. For a 3D images with a, b and c pixels in 3 dimensions, in `CACS.m` , the `Blocksize` is set as c, indicating number of 2D images;

2. To choose parameters for line-like or point-like signals, in `CalcLambda32.m` , keep the corresponding lambda formula, and comment another formula;

3. In `mainCS3D.m` , set XRes, YRes and ZRes as b, a and c
```
XRes = b;
YRes = a;
ZRes = c;
```
Note that if a and b are not equal (they are equal in the demo data) and their order are wrong, the program will not run correctly.

4. In `mainCS3D.m` , the enhance factor and voxel sizes for LR images are set in `PROCESSING PARAMETERS`
`FactorX/Y/Z` indicates the enhance factor and `xy/z sensor` indicates the voxel sizes for xy/z axis.

To run the CACS, open `CACS.m` and run.
The results will be saved at the same path as LR images.

To run the inference on your own data, make sure that:
1. Input the path of data in `CACS.m`.
2. The data can not be too large because it is limited by the memory of GPU, otherwise the data needs to be cropped with `crop_and_stitch.m`. With our GTX 1080 Ti, the largest 3D image has 100 pixels in three axis.

Then run the command as above. 
