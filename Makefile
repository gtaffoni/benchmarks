AUTH=gtaffoni
NAME=exabed-DGEMM
NNODES=2
CCODE=main
.DEFAULT_GOAL := help

help:
	@echo "Use \`make <target>\` where <target> is one of"
	@echo "  help     display this help message"
	@echo "  clean   clean the project"
	@echo "  run  CCODE=c_code  [NNODES=X]   execute the code c_code (defaul main.c) with optional numer of nodes (Default 2)"


#
OPT = -D_MMX_=8
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
VERBOSE      += -Wall -v
OBJECTS       = $(SOURCES:.c=.o)


#---------------------------------------------------------------------#
SYSTYPE      =  Exabed

#--------------------------- SYSTEM setup ----------------------------#


ifeq ($(SYSTYPE), 'exabed')
# OpenCL path


endif

ifeq ($(SYSTYPE), 'hotcat')
#
 CC        = mpicc
 CFLAGS    =
 LIBS      =
 LIBSDIR   =
 INCLUDE   =
 INCDIR    =
endif


.PHONY: all clean info test debug

all: $(PROG) $(OBJECTS) # $(ASSEMBLER)

$(PROG): $(OBJECTS)
	$(CC) $(DEBUG) $(OPTIMIZE) $^ -o $@  $(LIBS)
	@echo ' '
	@echo 'Program' $(PROG) 'compiled for' $(SYSTYPE) 'machine'
	@echo ' '

%.o: %.c
	$(CC) $(VERBOSE) $(DEBUG) $(OPTIMIZE) $(CFLAGS) $(OPT)  -I./include -c $< -o $@



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
