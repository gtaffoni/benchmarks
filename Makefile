AUTH=gtaffoni
NAME=exabed-DGEMM
NNODES=2
CCODE=main
.DEFAULT_GOAL := help

help:
	@echo "Use \`make <target>\` where <target> is one of"
	@echo "  help     display this help message"
	@echo "  clean   clean the project"
	@echo "  all     compile all .c files"
	@echo "  code_name compile codename.c"

#
OPT = -D_MMX_=4
#-----  Debug
OPT += -D_DEBUG_
#----------------------------------------------------------------------
#-----  OpenMP
OPT += -D_USE_OPENMP_ -fopenmp
#----------------------------------------------------------------------
#-----  OpenBlas
OPT += -D_OPEN_BLAS_

# Set the compiler options                                            #
#
OPTIMIZE      = -O3
OPTIMIZE     += -march=native # -Winline -finline-functions -Wextra
OPTIMIZE     += -ftree-vectorize
DEBUG        += # -ggdb
VERBOSE      += # -Wall -v
SOURCES       = $(wildcard *.c)
OBJECTS       = $(SOURCES:.c=.o)

#---------------------------------------------------------------------#
SYSTYPE      =  'gpu'

#--------------------------- SYSTEM setup ----------------------------#


ifeq ($(SYSTYPE), 'gpu')
CC        = nvc
CXX	  = nvc
CFLAGS    =
LIBS      =  -lgfortran
LIBSDIR   =
OPTIMIZE  = -O0
INCDIR    =
OPT       = -D_MMX_=4 -D_DEBUG_ -D_USE_OPENMP_ -mp
endif


#---------------------------------------------------------------------#
SYSTYPE      =  'hotcat'
#--------------------------- SYSTEM setup ----------------------------#



ifeq ($(SYSTYPE), 'hotcat')
CC        = mpicc
CXX				= gcc
CFLAGS    =
LIBS      = -lopenblas -lpthread -lgfortran
LIBSDIR   = -L$(OpenBLAS_LIB) -L$(MPI_LIB)
INCDIR    = -I$(MPI_INCLUDE) -I$(OpenBLAS_INC)
endif

ifeq ($(SYSTYPE), 'fpga')
CC        = mpiccc
CFLAGS    =
LIBS      =
LIBSDIR   = -L$(MPI_LIB)
INCDIR    = -I$(MPI_INCLUDE)
endif

.PHONY: all clean info test debug

PROG = dgemm_fpga dgemm_gpu dgemm_serial dgemm_open_blas

all: $(PROG) $(OBJECTS) # $(ASSEMBLER)


%.o: %.c
	$(CC) $(VERBOSE) $(DEBUG) $(OPTIMIZE) $(CFLAGS) $(OPT)  $(INCDIR)  -c $< -o $@


dgemm_open_blas: dgemm_open_blas.o Makefile
	$(CC) $(DEBUG) $(OPTIMIZE) $(LIBSDIR) $(OPT)  $< -o $@.x  $(LIBS)
	@echo 'Program compiled for' $(SYSTYPE) 'machine'

dgemm_open_blas_s: dgemm_open_blas_s.c Makefile
	$(CXX) $(VERBOSE) $(DEBUG) $(OPTIMIZE) $(CFLAGS) $(OPT)  $(INCDIR) $(LIBSDIR) $(LIBS)  dgemm_open_blas_s.c -o dgemm_open_blas_s.x
	@echo 'Program compiled for' $(SYSTYPE) 'machine'


dgemm_fpga: dgemm_fpga.o Makefile
	$(CC) $(DEBUG) $(OPTIMIZE) $(LIBSDIR) $(OPT)  $< -o $@.x  $(LIBS)
	@echo 'Program compiled for' $(SYSTYPE) 'machine'

dgemm_gpu: dgemm_gpu.c Makefile
		$(CXX) $(VERBOSE) $(DEBUG) $(OPTIMIZE) $(CFLAGS) $(OPT)  $(INCDIR) dgemm_gpu.c -o dgemm_gpu.x
		@echo 'Program compiled for' $(SYSTYPE) 'machine'

dgemm_serial: dgemm_serial.c Makefile
	$(CXX) $(VERBOSE) $(DEBUG) $(OPTIMIZE) $(CFLAGS) $(OPT)  $(INCDIR) dgemm_serial.c -o dgemm_serial.x
	@echo 'Program compiled for' $(SYSTYPE) 'machine'
clean:
	rm -rf $(PROG) $(OBJECTS) $(ASSEMBLER)

test: $(PROG)
	@echo 'oooOOO... testing ...OOOooo'
	mpirun  -np $(NNODES) ./$<
	@echo 'oooOOO... testing ...OOOooo'

debug: $(PROG)
	@echo 'oooOOO ... debugging ...OOOooo'
	gdb ./$< --args
	@echo 'oooOOO ... debugging ...OOOooo'
