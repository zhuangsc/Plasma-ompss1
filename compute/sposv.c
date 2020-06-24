/**
 *
 * @file sposv.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:45:07 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_sposv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N symmetric positive definite (or Hermitian positive definite
 *  in the complex case) matrix and X and B are N-by-NRHS matrices.
 *  The Cholesky decomposition is used to factor A as
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo = PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and  L is a lower triangular matrix.
 *  The factored form of A is then used to solve the system of equations A * X = B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_sposv_Tile
 * @sa PLASMA_sposv_Tile_Async
 * @sa PLASMA_cposv
 * @sa PLASMA_dposv
 * @sa PLASMA_sposv
 *
 ******************************************************************************/
int PLASMA_sposv(PLASMA_enum uplo, int N, int NRHS,
                 float *A, int LDA,
                 float *B, int LDB)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sposv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_sposv", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_sposv", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_sposv", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_sposv", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_sposv", "illegal value of LDB");
        return -7;
    }
    /* Quick return - currently NOT equivalent to LAPACK's
     * LAPACK does not have such check for DPOSV */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_SPOSV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sposv", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,    sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_sooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,   
                            sequence, &request);
        plasma_siplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_sposv_Tile_Async(uplo, &descA, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);
        plasma_sooptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);
        plasma_siptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        plasma_dynamic_sync();
    }
    
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_sposv_Tile - Solves a symmetric positive definite or Hermitian positive definite
 *  system of linear equations using the Cholesky factorization.
 *  Tile equivalent of PLASMA_sposv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_sposv
 * @sa PLASMA_sposv_Tile_Async
 * @sa PLASMA_cposv_Tile
 * @sa PLASMA_dposv_Tile
 * @sa PLASMA_sposv_Tile
*
 ******************************************************************************/
int PLASMA_sposv_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sposv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sposv_Tile_Async(uplo, A, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sposv_Tile_Async - Solves a symmetric positive definite or Hermitian
 *  positive definite system of linear equations using the Cholesky factorization.
 *  Non-blocking equivalent of PLASMA_sposv_Tile().
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
 * @sa PLASMA_sposv
 * @sa PLASMA_sposv_Tile
 * @sa PLASMA_cposv_Tile_Async
 * @sa PLASMA_dposv_Tile_Async
 * @sa PLASMA_sposv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sposv_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sposv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sposv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sposv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sposv_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sposv_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_sposv_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_sposv_Tile", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -1);
    }
    /* Quick return - currently NOT equivalent to LAPACK's
     * LAPACK does not have such check for DPOSV */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/
    plasma_parallel_call_4(plasma_pspotrf,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pstrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaTrans : PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        float, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pstrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaNoTrans : PlasmaTrans,
        PLASMA_enum, PlasmaNonUnit,
        float, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
