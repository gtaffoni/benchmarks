# benchmarks

This is a set of codes to estimate the peak performance of ExaNeSt Exascale
prototype.

The prototype is based on Arm CPUs and FPGAs as accelerators and a custom
network.

The code a a simplified version of DGEMM to use on FPGA. It is not optimized
on purpose.  

It is based on MPI OpenMP.

CPU version is based on OpenBLAS.


## HowTo

Modify the makefile for your compiler and libraries.
```
ifeq ($(SYSTYPE), 'yoursystem')
 CC       = mpicc
CFLAGS    =
LIBS      =
LIBSDIR   = -L$(MPI_LIB)
INCDIR    = -I$(MPI_INCLUDE)
endif
```

Compile the code:

```
make SYSTYPE=yoursystem dgemm_xxxx

```
