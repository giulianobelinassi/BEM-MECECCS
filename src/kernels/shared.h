#ifndef SHARED_H
#define SHARED_H

#include <thrust/complex.h>
#include "cublas_v2.h"

extern "C"{            

__constant__ FREAL delta[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                       
extern FREAL* device_cx;
extern FREAL* device_cy;                                                              
extern FREAL* device_cz;
extern FREAL* device_cxm;
extern FREAL* device_cym;                                                               
extern FREAL* device_czm;
extern FREAL* device_etas;
extern FREAL* device_gi;
extern FREAL* device_ome;
extern FREAL* device_hestdiag;
extern FREAL* device_gestdiag;
extern int*   device_cone;        
extern thrust::complex<FREAL>* device_zh;
extern thrust::complex<FREAL>* device_zdfi;
extern thrust::complex<FREAL>* device_zfi;
extern int*   device_kode;

extern int swapped;

void cuda_assert(cudaError_t error);

void cublas_assert(cublasStatus_t error);

int largest_possible_width(size_t sizeof_column_mem, int columns, int* iterations);
}

/*
https://stackoverflow.com/questions/24095248/cudas-warp-shuffle-for-double-precision-data
*/
__device__ __inline__ double shfl(double x, int lane)
{
    const int warpSize = 32;
    // Split the double number into 2 32b registers.
    int lo, hi;
    asm volatile("mov.b64 {%0,%1}, %2;":"=r"(lo),"=r"(hi):"d"(x));
    // Shuffle the two 32b registers.
    lo = __shfl_xor(lo,lane,warpSize);
    hi = __shfl_xor(hi,lane,warpSize);
    // Recreate the 64b number.
    asm volatile("mov.b64 %0,{%1,%2};":"=d"(x):"r"(lo),"r"(hi));
    return x;
}

__device__ __inline__ float shfl(float x, int lane)
{
	return __shfl_down(x, lane);
}

#endif

