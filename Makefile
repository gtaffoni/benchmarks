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
OPT = -D_MMX_=256
#-----  Debug
OPT += # -D_DEBUG_
#----------------------------------------------------------------------
#-----  OpenMP
OPT += -D_USE_OPENMP_ -fopenmp
#----------------------------------------------------------------------
#-----  OpenBlas
OPT += # -D_OPEN_BLAS_

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
SYSTYPE      =  'hotcat'

#--------------------------- SYSTEM setup ----------------------------#



ifeq ($(SYSTYPE), 'hotcat')
CC        = mpicc
CFLAGS    =
LIBS      = -lopenblas -lpthread -lgfortran
LIBSDIR   = -L$(OpenBLAS_LIB) -L$(MPI_LIB)
INCDIR    = -I$(MPI_INCLUDE) -I$(OpenBLAS_INC)
endif

ifeq ($(SYSTYPE), 'fpga')
CC        = mpicc
CFLAGS    =
LIBS      =
LIBSDIR   = -L$(MPI_LIB)
INCDIR    = -I$(MPI_INCLUDE)
endif

.PHONY: all clean info test debug

all: $(PROG) $(OBJECTS) # $(ASSEMBLER)

$(PROG): $(OBJECTS)
	$(CC) $(DEBUG) $(OPTIMIZE) $^ -o $@  $(LIBS)
	@echo ' '
	@echo 'Program' $(PROG) 'compiled for' $(SYSTYPE) 'machine'
	@echo ' '

%.o: %.c
	@echo $(SYSTYPE)
	$(CC) $(VERBOSE) $(DEBUG) $(OPTIMIZE) $(CFLAGS) $(OPT)  $(INCDIR)  -c $< -o $@



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
