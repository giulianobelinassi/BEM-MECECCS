#include "shared.h"
#include <cstdio>

extern "C"{

FREAL* device_cx;
FREAL* device_cy;
FREAL* device_cz;
FREAL* device_cxm;
FREAL* device_cym;
FREAL* device_czm;
FREAL* device_etas;
FREAL* device_gi;
FREAL* device_ome;
int*   device_cone;

/*Disponivel apos a execucao de Ghmatece*/
FREAL* device_hestdiag;
FREAL* device_gestdiag;


/*Inputecd*/
thrust::complex<FREAL>* device_zdfi;
int* device_kode;

/*Disponivel apos a execucao de Ghmatecd*/
thrust::complex<FREAL>* device_zh;
thrust::complex<FREAL>* device_zfi;

/*indica se precisou fazer swap ao monstar zh*/
int swapped = 0;

void cuda_assert(cudaError_t error)
{
    if (error != cudaSuccess)
    {   
        fputs(cudaGetErrorString(cudaGetLastError()), stderr);
        putc('\n', stderr);
        exit(1);
    }   
}

/*https://stackoverflow.com/questions/13041399/equivalent-of-cudageterrorstring-for-cublas*/
void cublas_assert(cublasStatus_t error)
{
	switch (error)
    {
        case CUBLAS_STATUS_NOT_INITIALIZED:
            fputs("CUBLAS_STATUS_NOT_INITIALIZED\n", stderr);
			exit(1);

        case CUBLAS_STATUS_ALLOC_FAILED:
            fputs("CUBLAS_STATUS_ALLOC_FAILED\n", stderr);
			exit(1);
        
		case CUBLAS_STATUS_INVALID_VALUE:
            fputs("CUBLAS_STATUS_INVALID_VALUE\n", stderr);
			exit(1);

        case CUBLAS_STATUS_ARCH_MISMATCH:
            fputs("CUBLAS_STATUS_ARCH_MISMATCH\n", stderr);
			exit(1);

        case CUBLAS_STATUS_MAPPING_ERROR:
            fputs("CUBLAS_STATUS_MAPPING_ERROR\n", stderr);
			exit(1);

        case CUBLAS_STATUS_EXECUTION_FAILED:
            fputs("CUBLAS_STATUS_EXECUTION_FAILED\n", stderr);
			exit(1);

        case CUBLAS_STATUS_INTERNAL_ERROR:
            fputs("CUBLAS_STATUS_INTERNAL_ERROR\n", stderr);
			exit(1);
		
		case CUBLAS_STATUS_SUCCESS:
			return;
    }
	fputs("CUBLAS: Erro desconhecido.\n", stderr);
	exit(1);

}

void send_shared_data_to_gpu_(
        FREAL cx[],
        FREAL cy[],
        FREAL cz[],
        FREAL cxm[],
        FREAL cym[],
        FREAL czm[],
        FREAL etas[],
        FREAL gi[],
        FREAL ome[],
        int cone[],
        int* np, 
        int* npg,
        int* n,
        int* nbe 
        )
{
    cudaError_t error;

    /*Aloque memória para os vetores na GPU*/
    error = cudaMalloc(&device_cone, 4*(*n)*sizeof(int));
    cuda_assert(error);

    error = cudaMalloc(&device_cx, (*np)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_cy, (*np)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_cz, (*np)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_cxm, (*n)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_cym, (*n)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_czm, (*n)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_gi, (*npg)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_ome, (*npg)*sizeof(FREAL));
    cuda_assert(error);

    error = cudaMalloc(&device_etas, (*n)*3*sizeof(FREAL));
    cuda_assert(error);

    /*mova os dados para lá*/
    error = cudaMemcpy(device_cone, cone, 4*(*n)*sizeof(int), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cx, cx, (*np)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cy, cy, (*np)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cz, cz, (*np)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cxm, cxm, (*n)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cym, cym, (*n)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_czm, czm, (*n)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_gi, gi, (*npg)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_ome, ome, (*npg)*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_etas, etas, (*n)*3*sizeof(FREAL), cudaMemcpyHostToDevice);
    cuda_assert(error);

}


void send_dyn_data_gpu_(int* nbe, thrust::complex<FREAL> zdfi[], int kode[])
{
    int nn = 3*(*nbe);
    cudaError_t error;

    error = cudaMalloc(&device_zdfi, nn*sizeof(thrust::complex<FREAL>));
    cuda_assert(error);
    error = cudaMalloc(&device_kode, nn*sizeof(int));
    cuda_assert(error);

    error = cudaMemcpy(device_zdfi, zdfi, nn*sizeof(thrust::complex<FREAL>), cudaMemcpyHostToDevice);
    cuda_assert(error);
    error = cudaMemcpy(device_kode, kode, nn*sizeof(int), cudaMemcpyHostToDevice);
    cuda_assert(error);
}

void deallocate_shared_gpu_data_()
{
    cudaError_t error;

    error = cudaFree(device_cone);
    cuda_assert(error);
    error = cudaFree(device_gi);
    cuda_assert(error);
    error = cudaFree(device_ome);
    cuda_assert(error);
    error = cudaFree(device_etas);
    cuda_assert(error);
    error = cudaFree(device_cx);
    cuda_assert(error);
    error = cudaFree(device_cz);
    cuda_assert(error);
    error = cudaFree(device_cxm);
    cuda_assert(error);
    error = cudaFree(device_cym);
    cuda_assert(error);
    error = cudaFree(device_czm);
    cuda_assert(error);
    error = cudaFree(device_hestdiag);
    cuda_assert(error);
    error = cudaFree(device_gestdiag);
    cuda_assert(error);
    error = cudaFree(device_zdfi);
    cuda_assert(error);
    error = cudaFree(device_kode);
    cuda_assert(error);
}

int largest_possible_width(size_t sizeof_column_mem, int columns, int* iterations)
{
	size_t available_mem;
	size_t total_mem;
	int possible_width;

	cuda_assert(cudaMemGetInfo(&available_mem, &total_mem));
//	available_mem = 8L*1024L*1024L*1024L; //Simulate a GPU with 8Gb of video memory

	if ((3*columns)*sizeof_column_mem < available_mem)
	{	*iterations = 1;
		return columns;
	}
	possible_width = available_mem/(3*sizeof_column_mem);

// As vezes o programa crasha pois, mesmo fazendo os cálculos exatos,
// A GPU fica sem memória. Sendo assim, eu retiro algumas colunas para
// tentar amenizar o problema
	possible_width -= 10;

	*iterations = (columns + possible_width - 1)/possible_width;

	return possible_width;
}

}
