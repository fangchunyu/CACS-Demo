/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "cuda.h"
#include "cufft.h"
#include <iostream>
#include <math.h>
#include <stdio.h>


//Constrains
static const int MAX_THREADS_CUDA = 1024; //adjust it for your GPU. This is correct for a 2.0 architecture
static const int MAX_BLOCKS_CUDA = 65535;
static const int dimsImage = 3;//so thing can be set at co0mpile time


void convolution3DfftCUDAFull(float* im, int* imDim, float* kernel, int* kernelDim, float * convResult);



/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
	float* im; int* imDim; float* kernel; int* kernelDim; float * convResult;

    /* Retrieve the input data */
	im = (float*)mxGetData(prhs[0]);
	imDim = (int*)mxGetData(prhs[1]);
	kernel = (float*)mxGetData(prhs[2]);
	kernelDim = (int*)mxGetData(prhs[3]);


	/* Create an mxArray for the output data */
	const mwSize dims[] = { imDim[0], imDim[1], imDim[2]};
	plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	//Create pointer for output data
	convResult = (float*)mxGetData(plhs[0]);

	convolution3DfftCUDAFull(im, imDim, kernel, kernelDim, convResult);	

}



__device__ static const float PI_2 = 6.28318530717958620f;
__device__ static const float PI_1 = 3.14159265358979310f;

////////////////////////////////////////////////////////////////////////////////
// Modulate Fourier image of padded data by Fourier image of padded kernel
// and normalize by FFT size
////////////////////////////////////////////////////////////////////////////////
//Adapted from CUDA SDK examples
__device__ void mulAndScale(cufftComplex& a, const cufftComplex& b, const float& c)
{
	cufftComplex t = { c * (a.x * b.x - a.y * b.y), c * (a.y * b.x + a.x * b.y) };
	a = t;
};

__global__ void __launch_bounds__(MAX_THREADS_CUDA)  modulateAndNormalize_kernel(cufftComplex *d_Dst, cufftComplex *d_Src, long long int dataSize, float c)
{
	long long int i = (long long int)blockDim.x * (long long int)blockIdx.x + (long long int)threadIdx.x;
	long long int offset = (long long int)blockDim.x * (long long int)gridDim.x;
	while (i < dataSize)
	{

		cufftComplex a = d_Src[i];
		cufftComplex b = d_Dst[i];

		mulAndScale(a, b, c);
		d_Dst[i] = a;

		i += offset;
	}
};

//we use nearest neighbor interpolation to access FFT coefficients in the kernel
__global__ void __launch_bounds__(MAX_THREADS_CUDA)  modulateAndNormalizeSubsampled_kernel(cufftComplex *d_Dst, cufftComplex *d_Src, int kernelDim_0, int kernelDim_1, int kernelDim_2, int imDim_0, int imDim_1, int imDim_2, long long int datasize, float c)
{

	float r_0 = ((float)kernelDim_0) / ((float)imDim_0); //ratio between image size and kernel size to calculate access
	float r_1 = ((float)kernelDim_1) / ((float)imDim_1);
	float r_2 = ((float)kernelDim_2) / ((float)imDim_2);

	long long int i = (long long int)blockDim.x * (long long int)blockIdx.x + (long long int)threadIdx.x;
	long long int offset = (long long int)blockDim.x * (long long int)gridDim.x;
	int k_0, k_1, k_2;
	int aux;
	float auxExp, auxSin, auxCos;
	while (i < datasize)
	{
		//for each dimension we need to access k_i*r_i  i=0, 1, 2
		aux = 1 + imDim_2 / 2;
		k_2 = i % aux;
		aux = (i - k_2) / aux;
		k_1 = aux % imDim_1;
		k_0 = (aux - k_1) / imDim_1;

		cufftComplex b = d_Dst[i];

		//apply shift in fourier domain since we did not apply fftshift to kernel (so we could use the trick of assuming the kernel is padded with zeros and then just subsample FFT)
		/* This is how we would do it in Matlab (linear phase change)
		auxExp = k_0 * r_0;
		auxExp += k_1 * r_1;
		auxExp += k_2 * r_2;
		auxExp *= PI_1;
		auxSin = sin(auxExp);
		auxCos = cos(auxExp);
		auxExp = b.x * auxCos - b.y * auxSin;

		b.y = b.x * auxSin + b.y * auxCos;
		b.x = auxExp;
		*/

		//add the ratio to each dimension and apply nearest neighbor interpolation
		//k_2 = min((int)(r_2*(float)k_2 + 0.5f),kernelDim_2-1);//the very end points need to be interpolated as "ceiling" instead of round or we can get oout of bounds access
		//k_1 = min((int)(r_1*(float)k_1 + 0.5f),kernelDim_1-1);
		//k_0 = min((int)(r_0*(float)k_0 + 0.5f),kernelDim_0-1);
		k_2 = ((int)(r_2*(float)k_2 + 0.5f)) % kernelDim_2;//the very end points need to be interpolated as "ceiling" instead of round or we can get oout of bounds access
		k_1 = ((int)(r_1*(float)k_1 + 0.5f)) % kernelDim_1;
		k_0 = ((int)(r_0*(float)k_0 + 0.5f)) % kernelDim_0;
		//calculate new coordinate relative to kernel size
		aux = 1 + kernelDim_2 / 2;
		cufftComplex a = d_Src[k_2 + aux *(k_1 + kernelDim_1 * k_0)];

		if ((k_0 + k_1 + k_2) % 2 == 1)//after much debugging it seems the phase shift is 0 or Pi (nothing in between). In Matlab is a nice linear change as programmed above
		{
			a.x = -a.x;
			a.y = -a.y;
		}
		mulAndScale(a, b, c);

		//__syncthreads();//this actually slows down the code by a lot (0.1 sec for 512x512x512)
		d_Dst[i] = a;

		i += offset;
	}
};

//WARNING: for cuFFT the fastest running index is z direction!!! so pos = z + imDim[2] * (y + imDim[1] * x)
__global__ void __launch_bounds__(MAX_THREADS_CUDA) fftShiftKernel(float* kernelCUDA, float* kernelPaddedCUDA, int kernelDim_0, int kernelDim_1, int kernelDim_2, int imDim_0, int imDim_1, int imDim_2)
{
	int kernelSize = kernelDim_0 * kernelDim_1 * kernelDim_2;

	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid<kernelSize)
	{
		//find coordinates
		long long int x, y, z, aux;


		z = tid % kernelDim_2;
		aux = (tid - z) / kernelDim_2;
		y = aux % kernelDim_1;
		x = (aux - y) / kernelDim_1;
		/*
		x = tid % kernelDim_0;
		aux = (tid - x)/kernelDim_0;
		y = aux % kernelDim_1;
		z = (aux - y)/kernelDim_1;
		*/
		//center coordinates
		x -= kernelDim_0 / 2;
		y -= kernelDim_1 / 2;
		z -= kernelDim_2 / 2;

		//circular shift if necessary
		if (x<0) x += imDim_0;
		if (y<0) y += imDim_1;
		if (z<0) z += imDim_2;

		//calculate position in padded kernel
		aux = z + imDim_2 * (y + imDim_1 * x);

		//aux = x + imDim_0 * (y + imDim_1 * z);

		//copy value
		kernelPaddedCUDA[aux] = kernelCUDA[tid];//for the most part it should be a coalescent access in oth places
	}
}


//=====================================================================
//WARNING: for cuFFT the fastest running index is z direction!!! so pos = z + imDim[2] * (y + imDim[1] * x)
//NOTE: to avoid transferring a large padded kernel, since memcpy is a limiting factor 
float* convolution3DfftCUDA(float* im, int* imDim, float* kernel, int* kernelDim)
{
	float* convResult = NULL;
	float* imCUDA = NULL;
	float* kernelCUDA = NULL;
	float* kernelPaddedCUDA = NULL;


	cufftHandle fftPlanFwd, fftPlanInv;


	long long int imSize = 1;
	long long int kernelSize = 1;
	for (int ii = 0; ii<dimsImage; ii++)
	{
		imSize *= (long long int) (imDim[ii]);
		kernelSize *= (long long int) (kernelDim[ii]);
	}

	long long int imSizeFFT = imSize + (long long int)(2 * imDim[0] * imDim[1]); //size of the R2C transform in cuFFTComplex

																				 //allocate memory for output result
	convResult = new float[imSize];

	//allocat ememory in GPU
	cudaMalloc((void**)&(imCUDA), imSizeFFT*sizeof(float));//a little bit larger to allow in-place FFT
	cudaMalloc((void**)&(kernelCUDA), (kernelSize)*sizeof(float));
	cudaMalloc((void**)&(kernelPaddedCUDA), imSizeFFT*sizeof(float));


	//TODO: pad image to a power of 2 size in all dimensions (use whatever  boundary conditions you want to apply)
	//TODO: pad kernel to image size
	//TODO: pad kernel and image to xy(z/2 + 1) for in-place transform
	//NOTE: in the example for 2D convolution using FFT in the Nvidia SDK they do the padding in the GPU, but in might be pushing the memory in the GPU for large images.

	//printf("Copying memory (kernel and image) to GPU\n");
	cudaMemcpy(kernelCUDA, kernel, kernelSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(imCUDA, im, imSize*sizeof(float), cudaMemcpyHostToDevice);

	
	//apply ffshift to kernel and pad it with zeros so we can calculate convolution with FFT
	cudaMemset(kernelPaddedCUDA, 0, imSizeFFT*sizeof(float));
	int numThreads = std::min((long long int)MAX_THREADS_CUDA, kernelSize);
	int numBlocks = std::min((long long int)MAX_BLOCKS_CUDA, (long long int)(kernelSize + (long long int)(numThreads - 1)) / ((long long int)numThreads));
	fftShiftKernel << <numBlocks, numThreads >> >(kernelCUDA, kernelPaddedCUDA, kernelDim[0], kernelDim[1], kernelDim[2], imDim[0], imDim[1], imDim[2]);


	//make sure GPU finishes before we launch two different streams
	cudaDeviceSynchronize();

	//printf("Creating R2C & C2R FFT plans for size %i x %i x %i\n",imDim[0],imDim[1],imDim[2]);
	cufftPlan3d(&fftPlanFwd, imDim[0], imDim[1], imDim[2], CUFFT_R2C);
	cufftSetCompatibilityMode(fftPlanFwd,CUFFT_COMPATIBILITY_NATIVE); //for highest performance since we do not need FFTW compatibility
	cufftPlan3d(&fftPlanInv, imDim[0], imDim[1], imDim[2], CUFFT_C2R);
	cufftSetCompatibilityMode(fftPlanInv,CUFFT_COMPATIBILITY_NATIVE);


	//transforming convolution kernel; TODO: if I do multiple convolutions with the same kernel I could reuse the results at teh expense of using out-of place memory (and then teh layout of the data is different!!!! so imCUDAfft should also be out of place)
	//NOTE: from CUFFT manual: If idata and odata are the same, this method does an in-place transform.
	//NOTE: from CUFFT manual: inplace output data xy(z/2 + 1) with fcomplex. Therefore, in order to perform an in-place FFT, the user has to pad the input array in the last dimension to Nn2 + 1 complex elements interleaved. Note that the real-to-complex transform is implicitly forward.
	cufftExecR2C(fftPlanFwd, imCUDA, (cufftComplex *)imCUDA);
	//transforming image
	cufftExecR2C(fftPlanFwd, kernelPaddedCUDA, (cufftComplex *)kernelPaddedCUDA);


	//multiply image and kernel in fourier space (and normalize)
	//NOTE: from CUFFT manual: CUFFT performs un-normalized FFTs; that is, performing a forward FFT on an input data set followed by an inverse FFT on the resulting set yields data that is equal to the input scaled by the number of elements.
	numThreads = std::min((long long int)MAX_THREADS_CUDA, imSizeFFT / 2);//we are using complex numbers
	numBlocks = std::min((long long int)MAX_BLOCKS_CUDA, (long long int)(imSizeFFT / 2 + (long long int)(numThreads - 1)) / ((long long int)numThreads));
	modulateAndNormalize_kernel << <numBlocks, numThreads >> >((cufftComplex *)imCUDA, (cufftComplex *)kernelPaddedCUDA, imSizeFFT / 2, 1.0f / (float)(imSize));//last parameter is the size of the FFT

																																								//inverse FFT 
	cufftExecC2R(fftPlanInv, (cufftComplex *)imCUDA, imCUDA);

	//copy result to host
	cudaMemcpy(convResult, imCUDA, sizeof(float)*imSize, cudaMemcpyDeviceToHost);

	//release memory
	(cufftDestroy(fftPlanInv));
	(cufftDestroy(fftPlanFwd));
	cudaFree(imCUDA);
	cudaFree(kernelCUDA);
	cudaFree(kernelPaddedCUDA);

	return convResult;
}


//WARNING: for cuFFT the fastest running index is z direction!!! so pos = z + imDim[2] * (y + imDim[1] * x)

void convolution3DfftCUDAFull(float* im, int* imDim, float* kernel, int* kernelDim, float * convResult) {

	int i, j, k, pos;
	float * kernelNew = (float *)malloc(kernelDim[0] * kernelDim[1] * kernelDim[2] * sizeof(float));

	//convert: send to convResult first
	for (k = 0; k<imDim[2]; k++) {
		for (j = 0; j<imDim[1]; j++) {
			for (i = 0; i<imDim[0]; i++) {
				pos = k + imDim[2] * (j + imDim[1] * i);
				convResult[pos] = im[k*imDim[0] * imDim[1] + j *imDim[0] + i];
			}
		}
	}

	//convert: kernel

	for (k = 0; k<kernelDim[2]; k++) {
		for (j = 0; j<kernelDim[1]; j++) {
			for (i = 0; i<kernelDim[0]; i++) {
				pos = k + kernelDim[2] * (j + kernelDim[1] * i);
				kernelNew[pos] = kernel[k*kernelDim[0] * kernelDim[1] + j *kernelDim[0] + i];			
			}
		}
	}
	

	//Do conv
	float* result = convolution3DfftCUDA(convResult, imDim, kernelNew, kernelDim);

	//Convert back
	for (k = 0; k<imDim[2]; k++) {
		for (j = 0; j<imDim[1]; j++) {
			for (i = 0; i<imDim[0]; i++) {
				pos = k + imDim[2] * (j + imDim[1] * i);
				convResult[k*imDim[0] * imDim[1] + j *imDim[0] + i] = result[pos];
			}
		}
	}


	free(kernelNew);
	free(result);

}

