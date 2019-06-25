#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>


#define TCPU_TIME (clock_gettime( id, &ts ), (double)ts.tv_sec +  \
                   (double)ts.tv_nsec * 1e-9)


void print_matrix(double **A, int *M, int *N);

int main(int argc, char *argv[]) {
        double **A, **B, **C;
        double alpha=1.0, beta=2.0;
        double *tmp;
        #ifdef _MMX_
        int NFPGA=_MMX_;
        #else
        int NFPGA=2;
        #endif

        struct timespec ts;
        clockid_t id = CLOCK_PROCESS_CPUTIME_ID;


        /* For OpenBlas compatibility */
        char ta = 'N';
        char tb = 'N';

        int i, j, k;
        int world_size=1;
        int N = world_size * NFPGA;
#ifdef _USE_OPENMP_
        int numThreads=1;
                #pragma omp parallel
        {
                #pragma omp Master
                numThreads = omp_get_num_threads();
        }
        printf ("Using %d OMP threads.\n", numThreads);

#endif
        /*
         * Memory allocation
         *
         * I want A, B, C to be contiguously allocated.
         * On slaves we can allocate less memory
         *
         */


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



        tmp = (double *) malloc (sizeof(double ) * N * N);
        B = (double **) malloc (sizeof(double *) * N);
        for (i = 0; i < N; i++)
                B[i] = &tmp[i * N];



        tmp = (double *) malloc (sizeof(double ) * N * N);
        C = (double **) malloc (sizeof(double *) * N);
        for (i = 0; i < N; i++)
                C[i] = &tmp[i * N];


        /*
         * Matrix initialization A=1., B=1.
         * could be done with random numbers
         */

#ifdef _USE_OPENMP_
                #pragma omp parallel for shared(A,B,C) private(i,j) schedule (dynamic, 10)
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





        /*
         * Begin Computation
         */

        printf("Begin DGEMM.\n");
        double startTime = TCPU_TIME;

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,NFPGA,N,N,alpha,*A, NFPGA, *B, N,beta,*C,N);

        double endTime = TCPU_TIME;

        printf("End Computation\n");
        double duration = endTime-startTime;
        /* VERIFY THIS!!!! */
        double gflops = 2.0 * N * N * N + 3.0 * N * N;
        printf("Execution Time =  %f\n", duration);
        printf("GFLOPs         =  %f\n", gflops/duration* 1.0e-9);



#ifdef _DEBUG_
        if ( N < 10) {
                print_matrix(C, &N, &N);
        }

#endif
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
