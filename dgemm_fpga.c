#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define T 11

int dgemm(int *M, int *N, int *K, double *alpha, double **A, double **B,
          double *beta, double **C, double **D);

void print_matrix(double **A, int *M, int *N);


/*
 *  Parallel DGEMM function.
 *
 *  D = alpha * A * B + beta * C
 *
 */

int main(int argc, char *argv[]) {
        double **A, **B, **C, **D;
        double alpha=1.0, beta=0.5;
        double *tmp;
        double startTime = 0.0;
        /* For OpenBlas compatibility */
        int world_rank, world_size;
        #ifdef _MMX_
        int NFPGA=_MMX_;
        #else
        int NFPGA=2;
        #endif
        int i;
        int j;

        MPI_Init(&argc, &argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        int N = world_size * NFPGA;
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
      #ifdef _USE_OPENMP_
                int numThreads=1;
                #pragma omp parallel
                {
                   #pragma omp Master
                   numThreads = omp_get_num_threads();
                }
                printf ("Using %d OMP threads.\n", numThreads);

        #endif
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
        B   = (double **) malloc (sizeof(double *) * N);
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


        if (world_rank == 0) {
                tmp = (double *) malloc (sizeof(double ) * N * N);
                D = (double **) malloc (sizeof(double *) * N);
                for (i = 0; i < N; i++)
                        D[i] = &tmp[i * N];
        }
        else {
                tmp = (double *) malloc (sizeof(double ) * N * NFPGA);
                D = (double **) malloc (sizeof(double *) * NFPGA);
                for (i = 0; i < NFPGA; i++)
                        D[i] = &tmp[i * N];
        }

        /*
         * Matrix initialization A=1., B=1. and C=1.
         * could be done with random numbers
         */
        if (world_rank == 0) {
#ifdef _USE_OPENMP_
        #pragma omp parallel for shared(A,B,C) private(i,j) schedule (static, 10)
#endif
                for (int i=0; i<N; i++) {
                        for (int j=0; j<N; j++) {
#ifdef _RANDOM_
                                A[i][j] = rand() % 101 - 50;;
                                B[i][j] = rand() % 101 - 50;;
                                C[i][j] = rand() % 101 - 50;;
#else
                                A[i][j] = 1.0;
                                B[i][j] = 1.0;
                                C[i][j] = 1.0;

#endif
                        }
                }
        }

        /* Initialize D --> D=0 everyone its own part*/
#ifdef _USE_OPENMP_
      #pragma omp parallel for shared(D) private(i,j) schedule (static, 10)
#endif
        for (i=0; i<NFPGA; i++) {
                for (j=0; j<N; j++) {
                        D[i][j] = 0.0;
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
                /* send A and C */
                int offset = NFPGA;
                int numElements = NFPGA * N;
                for (i=1; i<world_size; i++) {
                        MPI_Send(A[offset], numElements, MPI_DOUBLE, i, T, MPI_COMM_WORLD);
                        MPI_Send(C[offset], numElements, MPI_DOUBLE, i, T, MPI_COMM_WORLD);
                        offset += NFPGA;
                }
        }
        else {
                /* receive  A and C */
                MPI_Recv(A[0], NFPGA * N, MPI_DOUBLE, 0, T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(C[0], NFPGA * N, MPI_DOUBLE, 0, T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* Distribue B */
        MPI_Bcast(B[0], N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        dgemm(&NFPGA, &N, &N, &alpha, A, B, &beta, C, D);

/*  Master collects data from slaves, we are not using all_gather */
        if (world_rank == 0) {
                int offset = NFPGA;
                int numElements = NFPGA * N;
                for (i=1; i<world_size; i++) {
                        MPI_Recv(D[offset], numElements, MPI_DOUBLE, i, T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        offset += NFPGA;
                }
        }
        else {
                MPI_Send(D[0], NFPGA * N, MPI_DOUBLE, 0, T, MPI_COMM_WORLD);
        }

/* End of computation */
        if (world_rank == 0) {
                double endTime = MPI_Wtime();
                printf("End Computation\n");
                double duration = endTime-startTime;
                /* VERIFY THIS */
                /* DGEMM flop = 2*M*N*K + 3*M*N (see defintion of M, N, K below) in our case M=N=K*/
                double gflops = 2.0 * N * N * N + 3.0 * N * N;
                printf("Execution Time =  %f\n", duration);
                printf("GFLOPs         =  %f\n", gflops/duration * 1.0e-9);
        }


#ifdef _DEBUG_
        if (world_rank == 0 && N < 10) {
                print_matrix(D, &N, &N);
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
 *  DGEMM function.
 *
 *  D = alpha * A * B + beta * C
 *
 *
 *  Arguments
 *  ==========
 *
 *  M      - INTEGER.
 *           On entry,  M  specifies  the number  of rows  of the  matrix
 *            A,  of the  matrix  C and of the  matrix  D.
 *            M  must  be at least  2.
 *
 *
 *  N      - INTEGER.
 *           On entry,  N  specifies the number  of columns of the matrix
 *           B  and the number of columns of the matrix C. N must be
 *           at least 2.
 *
 *
 *  K      - INTEGER.
 *           On entry,  K  specifies  the number of columns of the matrix
 *           A  and the number of rows of the matrix B. K must
 *           be at least  2.
 *
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( M, K )
 *
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( K, N )
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *           supplied as zero then C need not be set on input.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( M, N ).
 *
 *  D      - DOUBLE PRECISION array of DIMENSION ( M, N ).
 *           On exit, the array  D  is overwritten by the  N by M  matrix
 *           ( alpha* A * B  + beta * C ).
 *
 */

int dgemm(int *M, int *N, int *K, double *alpha, double **A,
          double **B,  double *beta, double **C, double **D){

#ifdef _USE_OPENMP_
#pragma omp parallel for shared(A,B,C,D)  schedule (static, 10)
#endif
        for (int i=0; i<*M; i++) {
                for (int j=0; j<*N; j++) {
                        for (int kk=0; kk<*K; kk++) {
                                D[i][j] += A[i][kk] * B[kk][j];

                        }
                        D[i][j] = *alpha * D[i][j] + *beta * C[i][j];
                }
        }
        return(0);
}
