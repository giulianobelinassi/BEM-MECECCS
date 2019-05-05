# A Boundary Elements Method implementation for Stationary Elastodynamic Problems with GPU Acceleration

  Code implemented for a journal paper

## Quick Install:
Run the script
```
./deploy.sh
```

This will download and compile OpenBLAS and MAGMA. No sudo is required.
Inside the `src` directory, there is a script named `bench.sh` that reproduces the tests. The script `analyser.py` plot the time graphics.

### Dependencies

This project have mainly three dependencies: [OpenBLAS](http://www.openblas.net/), [CUDA](https://developer.nvidia.com/cuda-downloads) and [MAGMA](http://icl.cs.utk.edu/magma/). If you do not want to use GPU acceleration, then you only need OpenBLAS. Else all libraries are required.

#### 1. OpenBLAS

Download [OpenBLAS](http://www.openblas.net/) and untar it. Compile OpenBLAS with the following command:

```
INTERFACE64=1 USE_OPENMP=1 make -j <num_of_processors>
```
This will enable 64-bits indexing for very large matrices and OpenMP accelerated BLAS functions.

#### 2. MAGMA

Having OpenBLAS compiled, download [MAGMA](http://icl.cs.utk.edu/magma/) and unzip it. Although MAGMA supports a large set of BLAS implementations, we only tested our program with OpenBLAS.

Navigate to the unzipped MAGMA folder, then copy `make.inc-examples/make.inc.openblas` to `make.inc`. After that, find a block of text containing:

```
CFLAGS    = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -fopenmp
FFLAGS    = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument
F90FLAGS  = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input
NVCCFLAGS = -O3         -DNDEBUG -DADD_       -Xcompiler "$(FPIC)"
LDFLAGS   =     $(FPIC)                       -fopenmp
```
__Modify it to__
	
```
CFLAGS    = -O3 -DMAGMA_ILP64 $(FPIC) -DNDEBUG -DADD_ -Wall -fopenmp
FFLAGS    = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument
F90FLAGS  = -O3 $(FPIC) -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input
NVCCFLAGS = -O3 -DMAGMA_ILP64        -DNDEBUG -DADD_       -Xcompiler "$(FPIC)"
LDFLAGS   =     $(FPIC)                       -fopenmp
```

Also set the variable `OPENBLASDIR` with the path to the extracted and compiled OpenBLAS; and set `CUDADIR` to where you installed the CUDA library files. Adding `-DMAGMA_ILP64` enables 64-bit indexing for very large matrices. Then compile it with 
```
make dense
```

### Compiling the project

Currently we support both CPU-only and GPU-accelerated modes. The project only supports Linux AMD64 with __gcc >= 4.8__, __nvcc >= 8.0__ and __gfortran >= 4.8__.

#### 1. CPU-Only
For CPU-only, compile the OpenBLAS library as described above, update the variable `OPENBLASDIR` inside `Makefile` with the path to the compiled OpenBLAS folder, and compile the project with the command: `FRPEC=double make`. You can also compile with `FRPEC=float make` for single precision.

#### 2. GPU-accelerated

Compile both OpenBLAS and MAGMA as described above, update both variables `OPENBLASDIR` and `LIBMAGMADIR` inside `Makefile`wiith the path to the compiled OpenBLAS and MAGMA folders, respectively. Compile the project with the command `FRPEC=double make gpu`. You can also compile with `FRPEC=float make gpu` for single precision.

Currently, __We only support NVIDIA GPUs with Compute Capability > 3.0__. Check your card specifications. Note that some low-end GPUs are based on older NVIDIA architecture (like the GeForce GT 630M, that is Fermi-based), so make sure you have the required hardware.

### Executing the project.

```
./main <INPUT_STATIC> <INPUT_DYNAMIC> <OUTPUT_STATIC> <OUTPUT_DYNAMIC>
```
If the parameters provided are incorrect, the program will default to the files used while testing.

### Automated tests.
Run the script

```
./tests.sh
```
The tests requires the GPU-accelerated mode. 
### White paper
The paper is available on [ScienceDirect here](https://www.sciencedirect.com/science/article/abs/pii/S0955799718301395).
