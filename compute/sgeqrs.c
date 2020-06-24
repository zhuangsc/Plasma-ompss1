/**
 *
 * @file sgeqrs.c
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
 *  PLASMA_sgeqrs - Compute a minimum-norm solution min || A*X - B || using the
 *  QR factorization A = Q*R computed by PLASMA_sgeqrf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= M >= 0.
 *
 * @param[in] NRHS
 *          The number of columns of B. NRHS >= 0.
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned
 *          by PLASMA_sgeqrf.
 *          Remark: If out-of-place layout translation is used, the matrix A can
 *          be considered as input, however if inplace layout translation is
 *          enabled, the content of A will be reordered for computation and
 *          restored before exiting the function.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= M.
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by PLASMA_sgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgeqrs_Tile
 * @sa PLASMA_sgeqrs_Tile_Async
 * @sa PLASMA_cgeqrs
 * @sa PLASMA_dgeqrs
 * @sa PLASMA_sgeqrs
 * @sa PLASMA_sgeqrf
 *
 *******************************************************************************/
int PLASMA_sgeqrs(int M, int N, int NRHS,
                  float *A, int LDA,
                  PLASMA_desc *descT,
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
        plasma_fatal_error("PLASMA_sgeqrs", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_sgeqrs", "illegal value of M");
        return -1;
    }
    if (N < 0 || N > M) {
        plasma_error("PLASMA_sgeqrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_sgeqrs", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_sgeqrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, max(1, M))) {
        plasma_error("PLASMA_sgeqrs", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (min(M, min(N, NRHS)) == 0) {
        return PLASMA_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_SGELS, M, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeqrs", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N,    sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_sooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N,   
                            sequence, &request);
        plasma_siplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_sgeqrs_Tile_Async(&descA, descT, &descB, sequence, &request);

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
 *  PLASMA_sgeqrs_Tile - Computes a minimum-norm solution using the tile QR factorization.
 *  Tile equivalent of PLASMA_sgeqrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_sgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_sgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgeqrs
 * @sa PLASMA_sgeqrs_Tile_Async
 * @sa PLASMA_cgeqrs_Tile
 * @sa PLASMA_dgeqrs_Tile
 * @sa PLASMA_sgeqrs_Tile
 * @sa PLASMA_sgeqrf_Tile
 *
 ******************************************************************************/
int PLASMA_sgeqrs_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgeqrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sgeqrs_Tile_Async(A, T, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sgeqrs_Tile_Async - Computes a minimum-norm solution using the tile
 *  QR factorization.
 *  Non-blocking equivalent of PLASMA_sgeqrs_Tile().
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
 * @sa PLASMA_sgeqrs
 * @sa PLASMA_sgeqrs_Tile
 * @sa PLASMA_cgeqrs_Tile_Async
 * @sa PLASMA_dgeqrs_Tile_Async
 * @sa PLASMA_sgeqrs_Tile_Async
 * @sa PLASMA_sgeqrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sgeqrs_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgeqrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sgeqrs_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sgeqrs_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeqrs_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeqrs_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeqrs_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_sgeqrs_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, min(N, NRHS)) == 0) {
        return PLASMA_SUCCESS;
    }
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        plasma_parallel_call_7(plasma_psormqr,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaTrans,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_dynamic_call_8(plasma_psormqrrh,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaTrans,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    plasma_parallel_call_9(plasma_pstrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaUpper,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        float, 1.0,
        PLASMA_desc, plasma_desc_submatrix(descA, 0, 0, descA.n, descA.n),
        PLASMA_desc, plasma_desc_submatrix(descB, 0, 0, descA.n, descB.n),
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
