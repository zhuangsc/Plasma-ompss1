/**
 *
 * @file dsgesv.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated ds Tue Jan  7 11:45:09 2014
 *
 **/
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "common.h"

#define PLASMA_dlag2s(_descA, _descSB)                  \
  plasma_parallel_call_4(plasma_pdlag2s,                \
                         PLASMA_desc,      (_descA),    \
                         PLASMA_desc,      (_descSB),   \
                         PLASMA_sequence*, sequence,    \
                         PLASMA_request*,  request)

#define PLASMA_slag2d(_descSA, _descB)                  \
  plasma_parallel_call_4(plasma_pslag2d,                \
                         PLASMA_desc,      (_descSA),   \
                         PLASMA_desc,      (_descB),    \
                         PLASMA_sequence*, sequence,    \
                         PLASMA_request*,  request)

#define PLASMA_dlange(_norm, _descA, _result, _work)    \
  _result = 0;                                          \
  plasma_parallel_call_6(plasma_pdlange,                \
                         PLASMA_enum,      (_norm),     \
                         PLASMA_desc,      (_descA),    \
                         double*,          (_work),     \
                         double*,          &(_result),  \
                         PLASMA_sequence*, sequence,    \
                         PLASMA_request*,  request);

#define PLASMA_dlacpy(_descA, _descB)                           \
  plasma_parallel_call_5(plasma_pdlacpy,                        \
                         PLASMA_enum,      PlasmaUpperLower,    \
                         PLASMA_desc,      (_descA),            \
                         PLASMA_desc,      (_descB),            \
                         PLASMA_sequence*, sequence,            \
                         PLASMA_request*,  request)

#define PLASMA_dgeadd(_alpha, _descA, _descB)           \
  plasma_parallel_call_5(plasma_pdgeadd,                \
                         double, (_alpha),  \
                         PLASMA_desc,        (_descA),  \
                         PLASMA_desc,        (_descB),  \
                         PLASMA_sequence*,   sequence,  \
                         PLASMA_request*,    request)


#define PLASMA_dgetrs(_A, _IPIV, _X)                                    \
    {                                                                   \
        plasma_dynamic_call_3(                                          \
            plasma_pdbarrier_tl2pnl,                                    \
            PLASMA_desc, (_X),                                          \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*,  request);                                 \
                                                                        \
        plasma_dynamic_call_5(                                          \
            plasma_pdlaswp,                                             \
            PLASMA_desc, (_X),                                          \
            int *,       (_IPIV),                                       \
            int,         1,                                             \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*,  request);                                 \
                                                                        \
        plasma_parallel_call_9(                                         \
            plasma_pdtrsm,                                              \
            PLASMA_enum, PlasmaLeft,                                    \
            PLASMA_enum, PlasmaLower,                                   \
            PLASMA_enum, PlasmaNoTrans,                                 \
            PLASMA_enum, PlasmaUnit,                                    \
            double, 1.0,                                    \
            PLASMA_desc, (_A),                                          \
            PLASMA_desc, (_X),                                          \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*, request);                                  \
                                                                        \
        plasma_parallel_call_9(                                         \
            plasma_pdtrsm,                                              \
            PLASMA_enum, PlasmaLeft,                                    \
            PLASMA_enum, PlasmaUpper,                                   \
            PLASMA_enum, PlasmaNoTrans,                                 \
            PLASMA_enum, PlasmaNonUnit,                                 \
            double, 1.0,                                    \
            PLASMA_desc, (_A),                                          \
            PLASMA_desc, (_X),                                          \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*, request);                                  \
    }

#define PLASMA_sgetrs(_A, _IPIV, _X)                                    \
    {                                                                   \
        plasma_dynamic_call_3(                                          \
            plasma_psbarrier_tl2pnl,                                    \
            PLASMA_desc, (_X),                                          \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*,  request);                                 \
                                                                        \
        plasma_dynamic_call_5(                                          \
            plasma_pslaswp,                                             \
            PLASMA_desc, (_X),                                          \
            int *,       (_IPIV),                                       \
            int,         1,                                             \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*,  request);                                 \
                                                                        \
        plasma_parallel_call_9(                                         \
            plasma_pstrsm,                                              \
            PLASMA_enum, PlasmaLeft,                                    \
            PLASMA_enum, PlasmaLower,                                   \
            PLASMA_enum, PlasmaNoTrans,                                 \
            PLASMA_enum, PlasmaUnit,                                    \
            float, 1.0,                                    \
            PLASMA_desc, (_A),                                          \
            PLASMA_desc, (_X),                                          \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*, request);                                  \
                                                                        \
        plasma_parallel_call_9(                                         \
            plasma_pstrsm,                                              \
            PLASMA_enum, PlasmaLeft,                                    \
            PLASMA_enum, PlasmaUpper,                                   \
            PLASMA_enum, PlasmaNoTrans,                                 \
            PLASMA_enum, PlasmaNonUnit,                                 \
            float, 1.0,                                    \
            PLASMA_desc, (_A),                                          \
            PLASMA_desc, (_X),                                          \
            PLASMA_sequence*, sequence,                                 \
            PLASMA_request*, request);                                  \
    }

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dsgesv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 *
 *  PLASMA_dsgesv first attempts to factorize the matrix in COMPLEX and use this
 *  factorization within an iterative refinement procedure to produce a
 *  solution with COMPLEX*16 normwise backward error quality (see below).
 *  If the approach fails the method switches to a COMPLEX*16
 *  factorization and solve.
 *
 *  The iterative refinement is not going to be a winning strategy if
 *  the ratio COMPLEX performance over COMPLEX*16 performance is too
 *  small. A reasonable strategy should take the number of right-hand
 *  sides and the size of the matrix into account. This might be done
 *  with a call to ILAENV in the future. Up to now, we always try
 *  iterative refinement.
 *
 *  The iterative refinement process is stopped if ITER > ITERMAX or
 *  for all the RHS we have: RNRM < N*XNRM*ANRM*EPS*BWDMAX
 *  where:
 *
 *  - ITER is the number of the current iteration in the iterative refinement process
 *  - RNRM is the infinity-norm of the residual
 *  - XNRM is the infinity-norm of the solution
 *  - ANRM is the infinity-operator-norm of the matrix A
 *  - EPS is the machine epsilon returned by DLAMCH('Epsilon').
 *
 *  Actually, in its current state (PLASMA 2.1.0), the test is slightly relaxed.
 *
 *  The values ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 respectively.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The N-by-N coefficient matrix A. This matrix is not modified.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] IPIV
 *          On exit, the pivot indices that define the permutations.
 *
 * @param[in] B
 *          The N-by-NRHS matrix of right hand side matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[out] X
 *          If return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDX
 *          The leading dimension of the array B. LDX >= max(1,N).
 *
 * @param[out] ITER
 *          The number of the current iteration in the iterative refinement process
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsgesv_Tile
 * @sa PLASMA_dsgesv_Tile_Async
 * @sa PLASMA_dsgesv
 * @sa PLASMA_dgesv
 *
 ******************************************************************************/
int PLASMA_dsgesv(int N, int NRHS,
                  double *A, int LDA, int * IPIV,
                  double *B, int LDB,
                  double *X, int LDX, int *ITER)
{
    int NB;
    int status;
    PLASMA_desc  descA;
    PLASMA_desc  descB;
    PLASMA_desc  descX;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsgesv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    *ITER = 0;

    /* Check input arguments */
    if (N < 0) {
        plasma_error("PLASMA_dsgesv", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_dsgesv", "illegal value of NRHS");
        return -2;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_dsgesv", "illegal value of LDA");
        return -4;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_dsgesv", "illegal value of LDB");
        return -7;
    }
    if (LDX < max(1, N)) {
        plasma_error("PLASMA_dsgesv", "illegal value of LDX");
        return -9;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DSGESV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsgesv", "plasma_tune() failed");
        return status;
    }

    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    /* DOUBLE PRECISION INITIALIZATION */
    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,    sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)) );
        plasma_ddesc_alloc(  descX, NB, NB, N, NRHS, 0, 0, N, NRHS, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)); plasma_desc_mat_free(&(descX)) );
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,
                            sequence, &request);
        plasma_diplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS,
                            sequence, &request);

        descX = plasma_desc_init(
            PlasmaRealDouble, NB, NB, (NB*NB),
            LDX, NRHS, 0, 0, N, NRHS);
        descX.mat = X;
    }

    /* Call the native interface */
    status = PLASMA_dsgesv_Tile_Async(&descA, IPIV, &descB, &descX, ITER, sequence, &request);

    if (status == PLASMA_SUCCESS) {
        if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
            plasma_dooptile2lap( descX, X, NB, NB, LDX, NRHS,  sequence, &request);
            plasma_dynamic_sync();
            plasma_desc_mat_free(&descA);
            plasma_desc_mat_free(&descB);
            plasma_desc_mat_free(&descX);
        } else {
            plasma_diptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);
            plasma_diptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
            plasma_diptile2lap( descX, X, NB, NB, LDX, NRHS,  sequence, &request);
            plasma_dynamic_sync();
        }
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}


/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dsgesv_Tile - Solves a system of linear equations using the tile LU factorization
 *  and mixed-precision iterative refinement.
 *  Tile equivalent of PLASMA_dsgesv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          - If the iterative refinement converged, A is not modified;
 *          - otherwise, it fell back to double precision solution,
 *          and then A contains the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[out] IPIV
 *          On exit, the pivot indices that define the permutations (not equivalent to LAPACK).
 *
 * @param[in] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *
 * @param[out] X
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[out] ITER
 *          The number of the current iteration in the iterative refinement process
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsgesv
 * @sa PLASMA_dsgesv_Tile_Async
 * @sa PLASMA_dsgesv_Tile
 * @sa PLASMA_dgesv_Tile
 *
 ******************************************************************************/
int PLASMA_dsgesv_Tile(PLASMA_desc *A, int *IPIV,
                       PLASMA_desc *B, PLASMA_desc *X, int *ITER)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsgesv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    status = PLASMA_dsgesv_Tile_Async(A, IPIV, B, X, ITER, sequence, &request);
    if (status != PLASMA_SUCCESS)
        return status;
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dsgesv_Tile_Async - Solves a system of linear equations using
 *  the tile LU factorization and mixed-precision iterative refinement.
 *  Non-blocking equivalent of PLASMA_dsgesv_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsgesv
 * @sa PLASMA_dsgesv_Tile
 * @sa PLASMA_dsgesv_Tile_Async
 * @sa PLASMA_dgesv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dsgesv_Tile_Async(PLASMA_desc *A, int *IPIV,
                             PLASMA_desc *B, PLASMA_desc *X, int *ITER,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    int N, NB;
    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descX;
    PLASMA_desc descR, descSA, descSX;
    plasma_context_t *plasma;
    double *work;

    const int itermax = 30;
    const double bwdmax = 1.0;
    const double negone = -1.0;
    const double one = 1.0;
    int iiter;
    double Anorm, cte, eps, Rnorm, Xnorm;
    *ITER=0;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsgesv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dsgesv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dsgesv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsgesv_Tile", "invalid first descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsgesv_Tile", "invalid third descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descB = *B;
    }
    if (plasma_desc_check(X) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsgesv_Tile", "invalid fourth descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descX = *X;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb || descX.nb != descX.mb) {
        plasma_error("PLASMA_dsgesv_Tile", "only square tiles supported");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    /* Set N, NRHS, NT  */
    N     = descA.m;
    NB    = descA.nb;

    work = (double *)plasma_shared_alloc(plasma, PLASMA_SIZE, PlasmaRealDouble);
    if (work == NULL) {
        plasma_error("PLASMA_dsgesv", "plasma_shared_alloc() failed");
        plasma_shared_free(plasma, work);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    plasma_ddesc_alloc( descR,  NB, NB, descB.m, descB.n, 0, 0, descB.m, descB.n,
                        plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR) );
    plasma_sdesc_alloc( descSA, NB, NB, descA.m, descA.n, 0, 0, descA.m, descA.n,
                        plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR);
                        plasma_desc_mat_free(&descSA) );
    plasma_sdesc_alloc( descSX, NB, NB, descX.m, descX.n, 0, 0, descX.m, descX.n,
                        plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR);
                        plasma_desc_mat_free(&descSA); plasma_desc_mat_free(&descSX) );

    /* Compute some constants */
    PLASMA_dlange(PlasmaInfNorm, descA, Anorm, work);
    eps = LAPACKE_dlamch_work('e');

    /* Convert B from double precision to single precision and store
       the result in SX. */
    PLASMA_dlag2s(descB, descSX);
    if (sequence->status != PLASMA_SUCCESS)
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Convert A from double precision to single precision and store
       the result in SA. */
    PLASMA_dlag2s(descA, descSA);
    if (sequence->status != PLASMA_SUCCESS)
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Clear IPIV and Lbdl */
    plasma_memzero(IPIV, N, PlasmaInteger);

    /* Compute the LU factorization of SA */
    plasma_dynamic_call_3(
        plasma_psbarrier_tl2pnl,
        PLASMA_desc, descSA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    plasma_dynamic_call_4(plasma_psgetrf_rectil,
        PLASMA_desc, descSA,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Solve the system SA*SX = SB */
    PLASMA_sgetrs(descSA, IPIV, descSX);

    /* Convert SX back to double precision */
    PLASMA_slag2d(descSX, descX);

    /* Compute R = B - AX. */
    PLASMA_dlacpy(descB, descR);
    plasma_parallel_call_9(plasma_pdgemm,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNoTrans,
        double, negone,
        PLASMA_desc, descA,
        PLASMA_desc, descX,
        double, one,
        PLASMA_desc, descR,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Check whether the NRHS normwise backward error satisfies the
       stopping criterion. If yes return. Note that ITER=0 (already set). */
    PLASMA_dlange(PlasmaMaxNorm, descX, Xnorm, work);
    PLASMA_dlange(PlasmaMaxNorm, descR, Rnorm, work);

    /* Wait for the end of Anorm, Xnorm and Bnorm computations */
    plasma_dynamic_sync();

    cte = Anorm*eps*sqrt((double) N)*bwdmax;
    if (Rnorm < Xnorm * cte){
        /* The NRHS normwise backward errors satisfy the
           stopping criterion. We are good to exit. */
        plasma_desc_mat_free(&descSA);
        plasma_desc_mat_free(&descSX);
        plasma_desc_mat_free(&descR);
        plasma_shared_free(plasma, work);
        return PLASMA_SUCCESS;
    }

    /* Iterative refinement */
    for (iiter = 0; iiter < itermax; iiter++){

        /* Convert R from double precision to single precision
           and store the result in SX. */
        PLASMA_dlag2s(descR, descSX);

        /* Solve the system SA*SX = SB */
        PLASMA_sgetrs(descSA, IPIV, descSX);

        /* Convert SX back to double precision and update the current
           iterate. */
        PLASMA_slag2d(descSX, descR);
        PLASMA_dgeadd(one, descR, descX);

        /* Compute R = B - AX. */
        PLASMA_dlacpy(descB,descR);
        plasma_parallel_call_9(plasma_pdgemm,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaNoTrans,
            double, negone,
            PLASMA_desc, descA,
            PLASMA_desc, descX,
            double, one,
            PLASMA_desc, descR,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        /* Check whether the NRHS normwise backward errors satisfy the
           stopping criterion. If yes, set ITER=IITER>0 and return. */
        PLASMA_dlange(PlasmaInfNorm, descX, Xnorm, work);
        PLASMA_dlange(PlasmaInfNorm, descR, Rnorm, work);

        /* Wait for the end of Xnorm and Bnorm computations */
        plasma_dynamic_sync();

        if (Rnorm < Xnorm * cte){
            /* The NRHS normwise backward errors satisfy the
               stopping criterion. We are good to exit. */
            *ITER = iiter;

            plasma_desc_mat_free(&descSA);
            plasma_desc_mat_free(&descSX);
            plasma_desc_mat_free(&descR);
            plasma_shared_free(plasma, work);
            return PLASMA_SUCCESS;
        }
    }

    /* We have performed ITER=itermax iterations and never satisified
       the stopping criterion, set up the ITER flag accordingly and
       follow up on double precision routine. */
    *ITER = -itermax - 1;

    plasma_desc_mat_free(&descSA);
    plasma_desc_mat_free(&descSX);
    plasma_desc_mat_free(&descR);
    plasma_shared_free(plasma, work);

    /* Single-precision iterative refinement failed to converge to a
       satisfactory solution, so we resort to double precision. */

    /* Clear IPIV and Lbdl */
    plasma_memzero(IPIV, N, PlasmaInteger);

    plasma_dynamic_call_4(plasma_pdgetrf_rectil,
        PLASMA_desc, descA,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    PLASMA_dlacpy(descB, descX);

    PLASMA_dgetrs(descA, IPIV, descX);

    return PLASMA_SUCCESS;
}
