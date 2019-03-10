/// CUDA compileable code platform specific to Linux systems. Language standard: Host: C++11; Device: CUDA C. Requires GCC and NVCC linked to CUDA 8.0+ libraries.
/// This library is released for public use as material supplementary to PhD thesis: Development of Specialized Non-Linear Inversion Algorithms, Basis Functions, Eikonal Solvers, and Their Integration for Use in Joint Seismic and Gravitational Tomographic Inversion
/// Author: Zagid Abatchev - abatchev@ucla.edu ; PI: Paul Davis. Department of Earth, Planetary, and Space Sciences; University of California, Los Angeles; Date of publication: March 2019.///
/// This library is published as is and may have unresolved bugs. Use at own discretion. Users are encouraged to adapt some or all of this code without restriction. Refer to Github page for updated source code, located at: https://github.com/abatchev/ZTM .

///GFM_GPU.cu is a standalone gravity perturbation field forward model, implemented in CUDA C. Uses system shared virtual memory for fast file access, requires NVCC and CUDA 8.0+ to compile and properly run.///
///WARNING: This is an experimental implementation, and requires properly setting inversion specific grid size and decimation factors di, dj, dk, fh, fv inside main(). ///

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <fstream>
#include <iostream>

__global__ void GF(int bfr, int di, int dj, int dk, int fh, int fv, float *g, float *rho, float * s){                               //CUDA device kernel function GF() for forward modeling perturbation gravity fields as a function of perturbed density fields.
	int n=blockDim.x*blockIdx.x+threadIdx.x; float c=66.7f/(fh*fh*fv); int ii=di*fh, jj=dj*fh;                                      //Set raster index from CUDA kernel thread and block ID.
	float i=(n%ii)*1.f/fh, j=((n/ii)%jj)*1.f/fh, k=(n/(ii*jj))*1.f/fv;                                                              //Set gravity render field grid axial indices i,j,k.
	if(fabs(rho[n])>.001f){                                                                                                         //If absolute density perturbation is greater than .001cc, iterate over grid.
		for(int ni=bfr; ni<di-bfr; ni++){for(int nj=bfr; nj<dj-bfr; nj++){                                                          //Iterate through all internal points in grid at least "bfr" axial indices away from a boundary.
			int m=ni+nj*di; float i0=ni, j0=nj, k0=s[m];                                                                            //Compute flattened linear index and convert axial indices to floats.
			atomicAdd(&g[m], c*rho[n]*(k0-k)/powf((i0-i)*(i0-i)+(j0-j)*(j0-j)+(k0-k)*(k0-k),1.5));}}}}                              //Perform atomic addition of perturbation gravity due to the intergrated mass of a grid block approximated at far field as a point mass perturbation.

int main(){int di=61, dj=61, dk=35, fh=1, fv=2, bfr=10, ii=(di-1)*fh+1, jj=(dj-1)*fh+1, kk=(dk-1)*fv+1, N=ii*jj*kk, M=di*dj, lr,ls;  //Declare inversion specific array dimensions di,dj,dk; horizontal and vertical grid point decimation fh, fv; boundary buffer size bfr (where gravity is not rendered due to edge effects).
	float *h_r=(float *)malloc(4*N); FILE *R=fopen("/dev/shm/R.dat","r"); lr=fread(h_r,4,N,R); fclose(R);                            //Allocate host memory for imported density perturbation field h_r. Import perturbation density field from shared memory directory in R.dat
	float *h_s=(float *)malloc(4*M); FILE *S=fopen("/dev/shm/S.dat","r"); ls=fread(h_s,4,M,S); fclose(S); ls+=lr;                    //Allocate host memory for imported gravity data h_s. Import perturbation density field from shared memory directory in S.dat
	float *h_g=(float *)malloc(4*M); for(int n=0;n<M;n++){h_g[n]=0;}                                                                 //Allocate host memory for rendered gravity field h_g. Initialize with zeros.
	float *r; cudaMalloc((void **)&r,4*N); cudaMemcpy(r,h_r,4*N,cudaMemcpyHostToDevice);                                             //Allocate device memory for imported density perturbation field R.
	float *g; cudaMalloc((void **)&g,4*M); cudaMemcpy(g,h_g,4*M,cudaMemcpyHostToDevice);                                             //Allocate device memory for rendered gravity field g.
	float *s; cudaMalloc((void **)&s,4*M); cudaMemcpy(s,h_s,4*M,cudaMemcpyHostToDevice);                                             //Allocate device memory for imported gravity data S.
	GF<<<N/ii,ii>>>(bfr,di,dj,dk,fh,fv,g,r,s);                                                                                       //Execute cuda kernel on the GPU
	cudaMemcpy(h_g,g,4*M,cudaMemcpyDeviceToHost);                                                                                    //Copy rendered field at pointer g from device to host memory (as h_g)
	FILE *G=fopen("/dev/shm/G.dat","w"); fwrite(h_g,4,M,G); fclose(G);                                                               //Open gravity field render file G.dat in system shared memory, write rendered array, close file.
	cudaFree(r); cudaFree(g); cudaFree(s); free(h_r); free(h_g); free(h_s);}                                                         //Free host and device memory.
