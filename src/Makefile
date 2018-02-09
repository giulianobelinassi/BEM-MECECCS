FC      = gfortran
NVCC    = nvcc

OPENBLASDIR=/tmp/OpenBLAS-0.2.20
LIBMAGMADIR=/tmp/magma-2.3.0


ifeq ($(FPREC), double)
NVPREC = -DFREAL=double
FPREC  = -fdefault-real-8 -fdefault-double-8
else
NVPREC = -DFREAL=float
FPREC  = 
endif

FP_FLAGS=-fdefault-real-8 -fdefault-double-8
FCFLAGS = $(EXTRA) $(FPREC) -fPIC -g -cpp -flto -Wall -Wno-unused-dummy-argument #-march=native -Ofast -funroll-loops
PARALLEL= -fopenmp
NVFLAGS = $(EXTRA) $(NVPREC) -I$(LIBMAGMADIR)/include -DMAGMA_ILP64 -Xcompiler -fPIC -arch=sm_35 -use_fast_math -O3 -Xptxas --opt-level=3 -Xptxas --allow-expensive-optimizations=true 
FLFLAGS = $(OPENBLASDIR)/libopenblas.a -fopenmp -L$(OPENBLASDIR)
RM      = rm -f

all: main

cpu: main

gpu: MODE = -DUSE_GPU
gpu: FLFLAGS += $(LIBMAGMADIR)/lib/libmagma.a -lcudart -lcublas -lcusparse -lstdc++ -L/usr/local/cuda/lib64 -L$(LIBMAGMADIR)/lib
gpu: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Normvec.o Linsolve.o Main.o Ghmatecd_cu.o Ghmatece_cu.o Interec1_cu.o Linsolve_cu.o shared.o 
	$(FC) $(PARALLEL) $(FCFLAGS) -o main $^ $(FLFLAGS)

main: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Normvec.o Linsolve.o Main.o 
	$(FC) $(PARALLEL) $(FCFLAGS) -o $@ $^ $(FLFLAGS)

Main.o: Main.for Inputece.for Ghmatecd.for Ghmatece.for
	$(FC) $(FCFLAGS) $(PARALLEL) $(MODE) -c $<

Linsolve_cu.o: kernels/Linsolve_cu.cu Linsolve.for
	$(NVCC) $(NVFLAGS) -c $<

Ghmatecd_cu.o: kernels/Ghmatecd_cu.cu Sing_de.for Nonsingd.for
	$(NVCC) $(NVFLAGS) -c $<

Ghmatece_cu.o: kernels/Ghmatece_cu.cu kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

Interec1_cu.o: kernels/Interec1_cu.cu kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

shared.o: kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

Linsolve.o: Linsolve.for
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

Normvec.o: Normvec.for
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Inputece.o: Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputecd.o: Inputecd.for
	$(FC) $(FCFLAGS) -c $<

Ghmatece.o: Ghmatece.for Singge.o Nonsinge.o
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

Ghmatecd.o: Ghmatecd.for Sing_de.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

Interec.o: Interec.for Sigmaec.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

Outputec.o: Outputec.for
	$(FC) $(FCFLAGS) -c $<

Sigmaec.o: Sigmaec.for Gauleg.o
	$(FC) $(FCFLAGS) -c $<

Sing_de.o: Sing_de.for Gauleg.o
	$(FC) $(FCFLAGS) -c $<

Singge.o: Singge.for Nonsinge.o
	$(FC) $(FCFLAGS) -c $<

Solfund.o: Solfund.for
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Solfundif.o: Solfundif.for
	$(FC) $(FCFLAGS) -c $<

Solfune.o: Solfune.for
	$(FC) $(FCFLAGS) -c $<

Gauleg.o: Gauleg.for
	$(FC) $(FCFLAGS) -c $<

Nonsingd.o: Nonsingd.for
	$(FC) $(FCFLAGS) -c $<

Nonsinge.o: Nonsinge.for
	$(FC) $(FCFLAGS) -c $<

clean:
	$(RM) main *.o fort.*
