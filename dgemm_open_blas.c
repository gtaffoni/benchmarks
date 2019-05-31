#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#define T 11

int dgemm(char *ta, char *tb, int *M, int *N, int *K, double *alpha, double **A, int *LDA, double **B,
  int *LDB, double *beta, double **C, int *LDC, double **D);
void print_matrix(double **A, int *M, int *N);

int main(int argc, char *argv[]) {
        double **A, **B, **C;
        double alpha=1.0, beta=1.0;
        double *tmp;
        double startTime;
        #ifdef _MMX_
          int NFPGA=_MMX_;
        #else
           int NFPGA=2;
        #endif
        /* For OpenBlas compatibility */
        char ta = 'N';
        char tb = 'N';
        int world_rank, world_size;
        int  i, j, k;
        int chunkSize = 10;

        MPI_Init(&argc, &argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        int N = world_size * NFPGA;
        int numThreads = omp_get_num_threads();

        /*
         * Memory allocation
         *
         * I want A, B, C, D to be contiguously allocated.
         * On slaves we can allocate less memory
         *
         */

        if (world_rank == 0) {
                printf("----------------------------------------\n");
                printf("DGEMM implementation = a * A * B + b * C\n");
                printf("----------------------------------------\n");
                printf("N Matrix = %d \n", N);
                printf("N Local  = %d \n", NFPGA);
                printf ("Using %d OMP threads.\n", numThreads);
                tmp = (double *) malloc (sizeof(double ) * N * N);
                A = (double **) malloc (sizeof(double *) * N);
                for (int i = 0; i < N; i++)
                        A[i] = &tmp[i * N];
        }
        else {
                tmp = (double *) malloc (sizeof(double ) * N * NFPGA);
                A = (double **) malloc (sizeof(double *) * NFPGA);
                for (i = 0; i < NFPGA; i++)
                        A[i] = &tmp[i * N];
        }


        tmp = (double *) malloc (sizeof(double ) * N * N);
        B = (double **) malloc (sizeof(double *) * N);
        for (i = 0; i < N; i++)
                B[i] = &tmp[i * N];


        if (world_rank == 0) {
                tmp = (double *) malloc (sizeof(double ) * N * N);
                C = (double **) malloc (sizeof(double *) * N);
                for (i = 0; i < N; i++)
                        C[i] = &tmp[i * N];
        }
        else {
                tmp = (double *) malloc (sizeof(double ) * N * NFPGA);
                C = (double **) malloc (sizeof(double *) * NFPGA);
                for (i = 0; i < NFPGA; i++)
                        C[i] = &tmp[i * N];
        }

        /*
         * Matrix initialization A=1., B=1. and D=1.
         * could be done with random numbers
         */
        if (world_rank == 0) {

                for (i=0; i<N; i++) {
                        for (j=0; j<N; j++) {
                                A[i][j] =  1. ;//rand() % 101 - 50;;
                                B[i][j] = 1. ;//rand() % 101 - 50;;
                        }
                }
        }

        /* Initialize C --> C=0 everyone its own part*/
        for (i=0; i<NFPGA; i++) {
                for (j=0; j<N; j++) {
                        C[i][j] = 0.0;
                }
        }


        /*
         * Begin Computation
         */
        if (world_rank == 0) {
                printf("Begin DGEMM.\n");
                startTime = MPI_Wtime();
        }


        if (world_rank == 0) {
                /* send A and D */
                int offset = NFPGA;
                int numElements = NFPGA * N;
                for (i=1; i<world_size; i++) {
                        MPI_Send(A[offset], numElements, MPI_DOUBLE, i, T, MPI_COMM_WORLD);
                        MPI_Send(D[offset], numElements, MPI_DOUBLE, i, T, MPI_COMM_WORLD);
                        offset += NFPGA;
                }
        }
        else {
                /* receive  A and D */
                MPI_Recv(A[0], NFPGA * N, MPI_DOUBLE, 0, T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(D[0], NFPGA * N, MPI_DOUBLE, 0, T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* Distribue B */
        MPI_Bcast(B[0], N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,&NFPGA,&N,&N,&alpha,A, &NFPGA, B, &N,&N,C,&N);
        dgemm(&ta, &tb, &M, &N, &N, &alpha, A, &N, B, &N, &beta, D, &N, C);

/*  Master collects data from slaves, we are not using all_gather */
        if (world_rank == 0) {
                int offset = NFPGA;
                int numElements = NFPGA * N;
                for (i=1; i<world_size; i++) {
                        MPI_Recv(C[offset], numElements, MPI_DOUBLE, i, T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        offset += NFPGA;
                }
        }
        else {
                MPI_Send(C[0], NFPGA * N, MPI_DOUBLE, 0, T, MPI_COMM_WORLD);
        }

/* End of computation */
        if (world_rank == 0) {
                double endTime = MPI_Wtime();
                printf("End Computation\n");
                double duration = endTime-startTime;
                /* VERIFY THIS!!!! */
                double gflops = 2.0 * N * N * N + 3.0 * N * N;
                printf("Execution Time =  %f\n", endTime-startTime);
                printf("GFLOPs         =  %f\n", gflops/duration* 1.0e-6);
        }


#ifdef DEBUG
        if (world_rank == 0 && N < 10) {
                 print_matrix(C, &N, &N);
        }
#endif
        MPI_Finalize();
        return 0;
}
/*
 * print out matrix
 */
void print_matrix(double **A, int *M, int *N){
    for (int i=0; i<*N; i++) {
          for (int j=0; j<*M; j++) {
                  printf("%f ", A[i][j]);
          }
          printf("\n");
  }
}


/*
 *  Arguments from OpenBLAS for reference
 *  ==========
 *
 *  TRANA - CHARACTER*1.
 *           On entry, TRANA specifies the form of op( A ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANA = 'N' or 'n',  op( A ) = A.
 *
 *              TRANA = 'T' or 't',  op( A ) = A'.
 *
 *              TRANA = 'C' or 'c',  op( A ) = A'.
 *
 *           Unchanged on exit.
 *
 *  TRANB - CHARACTER*1.
 *           On entry, TRANB specifies the form of op( B ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANB = 'N' or 'n',  op( B ) = B.
 *
 *              TRANB = 'T' or 't',  op( B ) = B'.
 *
 *              TRANB = 'C' or 'c',  op( B ) = B'.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry,  M  specifies  the number  of rows  of the  matrix
 *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry,  N  specifies the number  of columns of the matrix
 *           op( B ) and the number of columns of the matrix C. N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry,  K  specifies  the number of columns of the matrix
 *           op( A ) and the number of rows of the matrix op( B ). K must
 *           be at least  zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *           k  when  TRANA = 'N' or 'n',  and is  m  otherwise.
 *           Before entry with  TRANA = 'N' or 'n',  the leading  m by k
 *           part of the array  A  must contain the matrix  A,  otherwise
 *           the leading  k by m  part of the array  A  must contain  the
 *           matrix A.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. When  TRANA = 'N' or 'n' then
 *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *           least  max( 1, k ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
 *           n  when  TRANB = 'N' or 'n',  and is  k  otherwise.
 *           Before entry with  TRANB = 'N' or 'n',  the leading  k by n
 *           part of the array  B  must contain the matrix  B,  otherwise
 *           the leading  n by k  part of the array  B  must contain  the
 *           matrix B.
 *           Unchanged on exit.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in the calling (sub) program. When  TRANB = 'N' or 'n' then
 *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
 *           least  max( 1, n ).
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *           supplied as zero then C need not be set on input.
 *           Unchanged on exit.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *           Before entry, the leading  m by n  part of the array  C must
 *           contain the matrix  C,  except when  beta  is zero, in which
 *           case C need not be set on entry.
 *           On exit, the array  C  is overwritten by the  m by n  matrix
 *           ( alpha*op( A )*op( B ) + beta*C ).
 *
 *  LDC    - INTEGER.
 *           On entry, LDC specifies the first dimension of C as declared
 *           in  the  calling  (sub)  program.   LDC  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *
 */
