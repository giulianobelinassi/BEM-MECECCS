#include "shared.h"
#include "../include/magma.h"

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <cstdio>

extern "C"{

void lu_assert(int status)
{
    if (status != 0)
    {
        fputs("ERRO: Matriz ZH singular\n", stderr);
        exit(1);
    }
}

static void reorganize_zh(int nn, thrust::complex<FREAL> zh[])
{
    cudaError_t error;
    size_t zh_size = nn*nn*sizeof(thrust::complex<FREAL>);

    error = cudaFree(device_zh);
    cuda_assert(error);

    error = cudaMalloc(&device_zh, zh_size);
    cuda_assert(error);

    error = cudaMemcpy(device_zh, zh, zh_size, cudaMemcpyHostToDevice);
    cuda_assert(error);
}

static int zh_fits_in_memory(int nn, thrust::complex<FREAL> zh[])
{
    cudaError_t error;
    size_t zh_size = nn*nn*sizeof(thrust::complex<FREAL>);
    size_t available_mem;
    size_t total_mem;

    error = cudaMemGetInfo(&available_mem, &total_mem);
    cuda_assert(error);

    if (available_mem < zh_size)
        return 0;
    return 1;

}

/*Fortran routine coded in Linsolve.for*/
void linsolve_cpu_(int* nn, thrust::complex<FREAL> zh[], thrust::complex<FREAL> zfi[]);

void cuda_linsolve_(
			int* nn,
			int* n,
            thrust::complex<FREAL> zh[],
			thrust::complex<FREAL> zfi[]
		)
{
    cudaError_t error;

	magma_int_t status = 0;

    thrust::complex<FREAL> one(1.0f, 0.0f);
    thrust::complex<FREAL> zero(0., 0.);
    
	magma_int_t* piv = (magma_int_t*) calloc(*nn, sizeof(magma_int_t));

	if (!piv)
	{
		fputs("Erro: Memória Insuficiente", stderr);
		exit(1);
	}

    magma_init();

//    cudaMemcpy(zh, device_zh, (*nn)*(*nn)*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);


    if (zh_fits_in_memory(*nn, zh))
    {
        if (swapped)
            reorganize_zh(*nn, zh);


        if (sizeof(FREAL) == 8)
            magma_zgesv_gpu(*nn, 1, (magmaDoubleComplex_ptr) device_zh, *nn, piv, (magmaDoubleComplex_ptr) device_zfi, *nn, &status);
        else
            magma_cgesv_gpu(*nn, 1, (magmaFloatComplex_ptr) device_zh, *nn, piv, (magmaFloatComplex_ptr) device_zfi, *nn, &status);
        
        lu_assert(status);

        error = cudaMemcpy(zfi, device_zfi, (*nn)*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);
        cuda_assert(error);

        error = cudaFree(device_zfi);
        cuda_assert(error);

        error = cudaFree(device_zh);
        cuda_assert(error);

    }
    else
    {
        error = cudaMemcpy(zfi, device_zfi, (*nn)*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);
        cuda_assert(error);

        error = cudaFree(device_zfi);
        cuda_assert(error);

        error = cudaFree(device_zh);
        cuda_assert(error);

        if (sizeof(FREAL) == 8)
            magma_zgesv(*nn, 1, (magmaDoubleComplex_ptr) zh, *nn, piv, (magmaDoubleComplex_ptr) zfi, *nn, &status);
        else
            magma_cgesv(*nn, 1, (magmaFloatComplex_ptr) zh, *nn, piv, (magmaFloatComplex_ptr) zfi, *nn, &status);
//        linsolve_cpu_(nn, zh, zfi);
      
        lu_assert(status);
    }

	magma_finalize();
    
	free(piv);
}
}
