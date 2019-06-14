#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int dgemm(int *M, int *N, int *K, double *alpha, double **A, double **B,
          double *beta, double **C, double **D);

void print_matrix(double **A, int *M, int *N);


/*
 *  Serial DGEMM function.
 *
 *  D = alpha * A * B + beta * C
 *
 */

int main(int argc, char *argv[]) {
        double **A, **B, **C, **D;
        double alpha=1.0, beta=1.0;
        double *tmp;
        double startTime = 0.0;
        #ifdef _MMX_
        int NFPGA=_MMX_;
        #else
        int NFPGA=2;
        #endif
        int i;
        int j;



        int world_size=1;

        int N = world_size * NFPGA;
        /*
         * Memory allocation
         *
         * I want A, B, C, D to be contiguously allocated.
         * On slaves we can allocate less memory
         *
         */


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
        {
                tmp = (double *) malloc (sizeof(double ) * N * NFPGA);
                A = (double **) malloc (sizeof(double *) * NFPGA);
                for (i = 0; i < NFPGA; i++)
                        A[i] = &tmp[i * N];
        }

        {
                tmp = (double *) malloc (sizeof(double ) * N * N);
                B   = (double **) malloc (sizeof(double *) * N);
                for (i = 0; i < N; i++)
                        B[i] = &tmp[i * N];
        }

        {
                tmp = (double *) malloc (sizeof(double ) * N * NFPGA);
                C = (double **) malloc (sizeof(double *) * NFPGA);
                for (i = 0; i < NFPGA; i++)
                        C[i] = &tmp[i * N];
        }

        {
                tmp = (double *) malloc (sizeof(double ) * N * NFPGA);
                D = (double **) malloc (sizeof(double *) * NFPGA);
                for (i = 0; i < NFPGA; i++)
                        D[i] = &tmp[i * N];
        }

        /*
         * Matrix initialization A=1., B=1. and D=1.
         * could be done with random numbers
         */
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


        /* Initialize D --> D=0 everyone its own part*/
#ifdef _USE_OPENMP_
      #pragma omp parallel for shared(D) private(i,j) schedule (static, 10)
#endif
        for (i=0; i<NFPGA; i++) {
                for (j=0; j<N; j++) {
                        D[i][j] = 0.0;
                }
        }


        dgemm(&NFPGA, &N, &N, &alpha, A, B, &beta, D, C);


        printf("End Computation\n");
//                double duration = endTime-startTime;
        /* VERIFY THIS */
        /* DGEMM flop = 2*M*N*K + 3*M*N (see defintion of M, N, K below) in our case M=N=K*/
//                double gflops = 2.0 * N * N * N + 3.0 * N * N;
//                printf("Execution Time =  %f\n", duration);
//                printf("GFLOPs         =  %f\n", gflops/duration * 1.0e-6);



#ifdef _DEBUG_

        print_matrix(C, &N, &N);

#endif
        return 0;
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

/*
 *  * print out matrix
 *   */
void print_matrix(double **A, int *M, int *N){
        for (int i=0; i<*N; i++) {
                for (int j=0; j<*M; j++) {
                        printf("%f ", A[i][j]);
                }
                printf("\n");
        }
}
