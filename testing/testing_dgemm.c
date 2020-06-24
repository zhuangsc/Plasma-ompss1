/**
 *
 * @file testing_dgemm.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:45:18 2014
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "testing_dmain.h"

#undef COMPLEX
#define REAL

static int check_solution(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                          double alpha, double *A, int LDA,
                          double *B, int LDB,
                          double beta, double *Cref, double *Cplasma, int LDC);

int testing_dgemm(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 8) {
        USAGE("GEMM", "alpha beta M N K LDA LDB LDC",
              "   - alpha  : alpha coefficient\n"
              "   - beta   : beta coefficient\n"
              "   - M      : number of rows of matrices A and C\n"
              "   - N      : number of columns of matrices B and C\n"
              "   - K      : number of columns of matrix A / number of rows of matrix B\n"
              "   - LDA    : leading dimension of matrix A\n"
              "   - LDB    : leading dimension of matrix B\n"
              "   - LDC    : leading dimension of matrix C\n");
        return -1;
    }

    PLASMA_Set(PLASMA_TILE_SIZE, 128);
    double alpha = (double) atol(argv[0]);
    double beta = (double) atol(argv[1]);
    int M     = atoi(argv[2]);
    int N     = atoi(argv[3]);
    int K     = atoi(argv[4]);
    int LDA   = atoi(argv[5]);
    int LDB   = atoi(argv[6]);
    int LDC   = atoi(argv[7]);

    double eps;
    int info_solution;
    int i, j, ta, tb;
    int LDAxK = LDA*max(M,K);
    int LDBxN = LDB*max(K,N);
    int LDCxN = LDC*N;

    double *A      = (double *)malloc(LDAxK*sizeof(double));
	#pragma omp register( [LDAxK]A )
    double *B      = (double *)malloc(LDBxN*sizeof(double));
	#pragma omp register( [LDBxN]B )
    double *C      = (double *)malloc(LDCxN*sizeof(double));
	#pragma omp register( [LDCxN]C )
    double *Cinit  = (double *)malloc(LDCxN*sizeof(double));
	#pragma omp register( [LDCxN]Cinit )
    double *Cfinal = (double *)malloc(LDCxN*sizeof(double));
	#pragma omp register( [LDCxN]Cfinal )

    /* Check if unable to allocate memory */
    if ((!A)||(!B)||(!Cinit)||(!Cfinal)){
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');
    printf("\n");
    printf("------ TESTS FOR PLASMA DGEMM ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING DGEMM
     */

    /* Initialize A, B, C */
    LAPACKE_dlarnv_work(IONE, ISEED, LDAxK, A);
    LAPACKE_dlarnv_work(IONE, ISEED, LDBxN, B);
    LAPACKE_dlarnv_work(IONE, ISEED, LDCxN, C);

#ifdef COMPLEX
    for (ta=0; ta<3; ta++) {
        for (tb=0; tb<3; tb++) {
#else
    for (ta=0; ta<2; ta++) {
        for (tb=0; tb<2; tb++) {
#endif
            for ( i = 0; i < M; i++)
                for (  j = 0; j < N; j++)
                    Cinit[LDC*j+i] = C[LDC*j+i];
            for ( i = 0; i < M; i++)
                for (  j = 0; j < N; j++)
                    Cfinal[LDC*j+i] = C[LDC*j+i];

            /* PLASMA DGEMM */
            PLASMA_dgemm(trans[ta], trans[tb], M, N, K, alpha, A, LDA, B, LDB, beta, Cfinal, LDC);

            /* Check the solution */
            info_solution = check_solution(trans[ta], trans[tb], M, N, K, 
                                           alpha, A, LDA, B, LDB, beta, Cinit, Cfinal, LDC);

            if (info_solution == 0) {
                printf("***************************************************\n");
                printf(" ---- TESTING DGEMM (%s, %s) ............... PASSED !\n", transstr[ta], transstr[tb]);
                printf("***************************************************\n");
            }
            else {
                printf("************************************************\n");
                printf(" - TESTING DGEMM (%s, %s) ... FAILED !\n", transstr[ta], transstr[tb]);
                printf("************************************************\n");
            }
        }
    }
#ifdef _UNUSED_
    }}
#endif
    free(A); free(B); free(C);
    free(Cinit); free(Cfinal);

    return 0;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_solution(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                          double alpha, double *A, int LDA,
                          double *B, int LDB,
                          double beta, double *Cref, double *Cplasma, int LDC)
{
    int info_solution;
    double Anorm, Bnorm, Cinitnorm, Cplasmanorm, Clapacknorm, Rnorm, result;
    double eps;
    double beta_const;

    double *work = (double *)malloc(max(K,max(M, N))* sizeof(double));
    int Am, An, Bm, Bn;

    beta_const  = -1.0;

    if (transA == PlasmaNoTrans) {
        Am = M; An = K;
    } else {
        Am = K; An = M;
    }
    if (transB == PlasmaNoTrans) {
        Bm = K; Bn = N;
    } else {
        Bm = N; Bn = K;
    }

    Anorm       = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), Am, An, A,       LDA, work);
    Bnorm       = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), Bm, Bn, B,       LDB, work);
    Cinitnorm   = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M,  N,  Cref,    LDC, work);
    Cplasmanorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M,  N,  Cplasma, LDC, work);

    cblas_dgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, 
                (alpha), A, LDA, B, LDB, (beta), Cref, LDC);

    Clapacknorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Cref, LDC, work);

    cblas_daxpy(LDC * N, (beta_const), Cplasma, 1, Cref, 1);

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Cref, LDC, work);

    eps = LAPACKE_dlamch_work('e');

    printf("Rnorm %e, Anorm %e, Bnorm %e, Cinitnorm %e, Cplasmanorm %e, Clapacknorm %e\n", 
           Rnorm, Anorm, Bnorm, Cinitnorm, Cplasmanorm, Clapacknorm);

    result = Rnorm / ((Anorm + Bnorm + Cinitnorm) * N * eps);
    printf("============\n");
    printf("Checking the norm of the difference against reference DGEMM \n");
    printf("-- ||Cplasma - Clapack||_oo/((||A||_oo+||B||_oo+||C||_oo).N.eps) = %e \n", 
           result);

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
         printf("-- The solution is suspicious ! \n");
         info_solution = 1;
    }
    else {
         printf("-- The solution is CORRECT ! \n");
         info_solution= 0 ;
    }

    free(work);

    return info_solution;
}

int timing_dgemm(int argc, char **argv)
{

	int transa; int transb;
	int m, n, k;
	double alpha, beta;
	int lda, ldb, ldc;
	double *A, *B, *C;
	double *C0;
	int bs, rep;
	int i;

	double start, end;
	double elapsed;
	FILE *log;
	int num_threads;


	//gemm transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc
	if (argc != 12) {
		fprintf(stderr, "GEMMs transa transb m n k alpha lda ldb beta ldc bs rep\n");
		return 1;;
	}

	sscanf(argv[0], "%d", &transa);
	sscanf(argv[1], "%d", &transb);

	if ( !transa )
		transa = PlasmaNoTrans;
	else
		transa = PlasmaTrans;
	if ( !transb )
		transb = PlasmaNoTrans;
	else
		transb = PlasmaTrans;

	sscanf(argv[2], "%d", &m);
	sscanf(argv[3], "%d", &n);
	sscanf(argv[4], "%d", &k);
	sscanf(argv[5], "%lf", &alpha);
	sscanf(argv[6], "%d", &lda);
	sscanf(argv[7], "%d", &ldb);
	sscanf(argv[8], "%lf", &beta);
	sscanf(argv[9], "%d", &ldc);
	sscanf(argv[10], "%d", &bs);
	sscanf(argv[11], "%d", &rep);

	int dimA = max(m,k) * lda;
	int dimB = max(n,k) * ldb;
	int dimC = max(m,n) * ldc;
	A = malloc( dimA * sizeof(double));
	#pragma omp register ([dimA]A)
	B = malloc( dimB * sizeof(double));
	#pragma omp register ([dimB]B)
	C = malloc( dimC * sizeof(double));
	#pragma omp register ([dimC]C)

	LAPACKE_dlarnv(IONE, ISEED, dimA, A);
	LAPACKE_dlarnv(IONE, ISEED, dimB, B);
	LAPACKE_dlarnv(IONE, ISEED, dimC, C);

	PLASMA_Set( PLASMA_TILE_SIZE, bs);

	elapsed = 0.0;
	for ( i = 0; i < rep; i++ ) {
		start = gtime();
		PLASMA_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
		end = gtime();
		elapsed += end - start;
	}

	num_threads = omp_get_max_threads();
	dump_info("plasma_dgemm.log", num_threads, elapsed, rep);

	free(A); free(B); free(C);
    return 0;
}
