#include "shared.h"

#include <cstdio>
#include <cmath>
#include <thrust/complex.h>
#include <cuda_runtime.h>
//#include "cublas_v2.h"

extern "C"{

/*Reaproveite este kernel pois ele pode ser usado para resolver a parte
 * de deslocamentos de pontos internos de Interec.for*/
__global__ void ghmatecd_kernel(
                           int cone[],
                           FREAL cx[],
                           FREAL cy[],
                           FREAL cz[],
                           FREAL cxm[],
                           FREAL cym[],
                           FREAL czm[],
                           thrust::complex<FREAL> zh[],
                           thrust::complex<FREAL> zg[],
                           FREAL rn[][3],
                           thrust::complex<FREAL> zge,
                           thrust::complex<FREAL> zcs,
                           thrust::complex<FREAL> zcp,
                           FREAL fr,
                           FREAL gi[],
                           FREAL ome[],
                           FREAL c1,
                           FREAL c2,
                           FREAL c3,
                           FREAL c4,
						   FREAL hestdiag[],
						   FREAL gestdiag[],
                           int n,
                           int nbe,
                           int dim_cone,
						   int column_pad,
                           int fast_singular,
                           int* ret
                           );

void cuda_interec1_(int* n,
                    int* nbe,
                    int* npg,
                    int* l,
                    int* np,
					FREAL cxi[],
					FREAL cyi[],
					FREAL czi[],
                    thrust::complex<FREAL>* zge,
                    thrust::complex<FREAL>* zcs,
                    thrust::complex<FREAL>* zcp,
                    FREAL* c1,
                    FREAL* c2,
                    FREAL* c3,
                    FREAL* c4,
                    FREAL* fr,
                    thrust::complex<FREAL> zdfi[],
					thrust::complex<FREAL> zfi[],
                    thrust::complex<FREAL> zdsol[],
                    int* status
                   )
{
	
	dim3 threadsPerBlock(*npg,*npg);
	dim3 numBlocks(*nbe, *l);
	int shared_mem_size = 2*3*3*(*npg)*(*npg)*sizeof(thrust::complex<FREAL>);
	cudaError_t error;
    
	thrust::complex<FREAL> zhelem[3][3];
	thrust::complex<FREAL> zgelem[3][3];

	thrust::complex<FREAL>* device_zh;
	thrust::complex<FREAL>* device_zg;
    thrust::complex<FREAL>* device_zdfi;
    thrust::complex<FREAL>* device_zfi;
    thrust::complex<FREAL>* device_zdsol;

	FREAL* device_cxi;
	FREAL* device_cyi;
	FREAL* device_czi;

	int* device_return_status;
	int return_status;

    thrust::complex<FREAL> one(1., 0.);
    thrust::complex<FREAL> zero(0., 0.);
    thrust::complex<FREAL> minus_one(-1., 0.);

    cublasHandle_t handle;
	cublasStatus_t stats;

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	error = cudaMalloc(&device_zh, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	error = cudaMalloc(&device_zg, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	error = cudaMemset(device_return_status, 0, sizeof(int));
	cuda_assert(error);

	error = cudaMemset(device_zh, 0, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	error = cudaMemset(device_zg, 0, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	error = cudaMalloc(&device_cxi, (*l)*sizeof(FREAL));
	cuda_assert(error);
	
	error = cudaMalloc(&device_cyi, (*l)*sizeof(FREAL));
	cuda_assert(error);
	
	error = cudaMalloc(&device_czi, (*l)*sizeof(FREAL));
	cuda_assert(error);

	error = cudaMemcpy(device_cxi, cxi, (*l)*sizeof(FREAL), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cyi, cyi, (*l)*sizeof(FREAL), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_czi, czi, (*l)*sizeof(FREAL), cudaMemcpyHostToDevice);
	cuda_assert(error);

	ghmatecd_kernel<<<numBlocks, threadsPerBlock, shared_mem_size>>>(
						device_cone,
						device_cx,
						device_cy,
						device_cz,
						device_cxi,
						device_cyi,
						device_czi,
                        device_zh,
                        device_zg,
						(FREAL (*)[3]) device_etas,
						*zge,
						*zcs,
						*zcp,
						*fr,
						device_gi,
						device_ome,
						*c1,
						*c2,
						*c3,
						*c4,
						NULL,
						NULL,
						*nbe,
						*l,
                        *n,
						0,
						0,
						device_return_status
						);
    
	cudaDeviceSynchronize();

	error = cudaMalloc(&device_zdsol, 3*(*l)*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

    error = cudaMalloc(&device_zdfi, 3*(*nbe)*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);
	
	error = cudaMalloc(&device_zfi, 3*(*nbe)*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	error = cudaMemcpy(device_zdfi, zdfi, 3*(*nbe)*sizeof(thrust::complex<FREAL>), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_zfi, zfi, 3*(*nbe)*sizeof(thrust::complex<FREAL>), cudaMemcpyHostToDevice);
	cuda_assert(error);
    cudaDeviceSynchronize();

	error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	if (return_status != 0)
	{
		fputs("Matriz Singular\n", stderr);
	}

    stats = cublasCreate(&handle);
	cublas_assert(stats);
  
    if (sizeof(FREAL) == 8)
    {
        stats = cublasZgemv(handle, CUBLAS_OP_N, 3*(*l), 3*(*nbe), (cuDoubleComplex*) &one, (cuDoubleComplex*) device_zg, 3*(*l), (cuDoubleComplex*) device_zdfi, 1, (cuDoubleComplex*) &zero, (cuDoubleComplex*) device_zdsol, 1);
        cublas_assert(stats);
        cudaDeviceSynchronize();
        
        stats = cublasZgemv(handle, CUBLAS_OP_N, 3*(*l), 3*(*nbe), (cuDoubleComplex*) &(minus_one), (cuDoubleComplex*) device_zh, 3*(*l), (cuDoubleComplex*) device_zfi, 1, (cuDoubleComplex*) &one, (cuDoubleComplex*) device_zdsol, 1);
        cublas_assert(stats);
        cudaDeviceSynchronize(); 
    }
    else
    {
        stats = cublasCgemv(handle, CUBLAS_OP_N, 3*(*l), 3*(*nbe), (cuComplex*) &one, (cuComplex*) device_zg, 3*(*l), (cuComplex*) device_zdfi, 1, (cuComplex*) &zero, (cuComplex*) device_zdsol, 1);
        cublas_assert(stats);
        cudaDeviceSynchronize();
        
        stats = cublasCgemv(handle, CUBLAS_OP_N, 3*(*l), 3*(*nbe), (cuComplex*) &(minus_one), (cuComplex*) device_zh, 3*(*l), (cuComplex*) device_zfi, 1, (cuComplex*) &one, (cuComplex*) device_zdsol, 1);
        cublas_assert(stats);
        cudaDeviceSynchronize(); 
    }

	error = cudaMemcpy(zdsol, device_zdsol, 3*(*l)*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);
	cuda_assert(error);

    error = cudaFree(device_zh);
	cuda_assert(error);
	error = cudaFree(device_zg);
	*status = return_status;
	error = cudaFree(device_return_status);
	cuda_assert(error);

	error = cudaFree(device_cxi);
	cuda_assert(error);
	error = cudaFree(device_cyi);
	cuda_assert(error);
	error = cudaFree(device_czi);
	cuda_assert(error);
	error = cudaFree(device_zfi);
	cuda_assert(error);
	error = cudaFree(device_zdfi);
	cuda_assert(error);

	error = cudaFree(device_zdsol);
	cuda_assert(error);

	cublasDestroy(handle);
}
}
