/**
 *
 * @file sgetrs_incpiv.c
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
 *  PLASMA_sgetrs_incpiv - Solves a system of linear equations A * X = B, with a general N-by-N matrix A
 *  using the tile LU factorization computed by PLASMA_sgetrf_incpiv.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = PlasmaNoTrans:   A * X = B     (No transpose)
 *          = PlasmaTrans:     A**T * X = B  (Transpose)
 *          = PlasmaTrans: A**T * X = B  (Conjugate transpose)
 *          Currently only PlasmaNoTrans is supported.
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by PLASMA_sgetrf_incpiv.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] descL
 *          Auxiliary factorization data, related to the tile L factor, computed by PLASMA_sgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from PLASMA_sgetrf_incpiv (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgetrs_incpiv_Tile
 * @sa PLASMA_sgetrs_incpiv_Tile_Async
 * @sa PLASMA_cgetrs_incpiv
 * @sa PLASMA_dgetrs_incpiv
 * @sa PLASMA_sgetrs_incpiv
 * @sa PLASMA_sgetrf_incpiv
 *
 ******************************************************************************/
int PLASMA_sgetrs_incpiv(PLASMA_enum trans, int N, int NRHS,
                         float *A, int LDA,
                         PLASMA_desc *descL, const int *IPIV,
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
        plasma_fatal_error("PLASMA_sgetrs_incpiv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (trans != PlasmaNoTrans) {
        plasma_error("PLASMA_sgetrs_incpiv", "only PlasmaNoTrans supported");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (N < 0) {
        plasma_error("PLASMA_sgetrs_incpiv", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_sgetrs_incpiv", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_sgetrs_incpiv", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_sgetrs_incpiv", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_SGESV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgetrs_incpiv", "plasma_tune() failed");
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
    PLASMA_sgetrs_incpiv_Tile_Async(&descA, descL, IPIV, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
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
 *  PLASMA_sgetrs_incpiv_Tile - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Tile equivalent of PLASMA_sgetrs_incpiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by PLASMA_sgetrf_incpiv.
 *
 * @param[in] L
 *          Auxiliary factorization data, related to the tile L factor, computed by PLASMA_sgetrf_incpiv.
 *
 * @param[in] IPIV
 *          The pivot indices from PLASMA_sgetrf_incpiv (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgetrs_incpiv
 * @sa PLASMA_sgetrs_incpiv_Tile_Async
 * @sa PLASMA_cgetrs_incpiv_Tile
 * @sa PLASMA_dgetrs_incpiv_Tile
 * @sa PLASMA_sgetrs_incpiv_Tile
 * @sa PLASMA_sgetrf_incpiv_Tile
 *
 ******************************************************************************/
int PLASMA_sgetrs_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, const int *IPIV, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgetrs_incpiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sgetrs_incpiv_Tile_Async(A, L, IPIV, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sgetrs_incpiv_Tile_Async - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Non-blocking equivalent of PLASMA_sgetrs_incpiv_Tile().
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
 * @sa PLASMA_sgetrs_incpiv
 * @sa PLASMA_sgetrs_incpiv_Tile
 * @sa PLASMA_cgetrs_incpiv_Tile_Async
 * @sa PLASMA_dgetrs_incpiv_Tile_Async
 * @sa PLASMA_sgetrs_incpiv_Tile_Async
 * @sa PLASMA_sgetrf_incpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sgetrs_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, const int *IPIV, PLASMA_desc *B,
                                    PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descL;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgetrs_incpiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sgetrs_incpiv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sgetrs_incpiv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgetrs_incpiv_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(L) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgetrs_incpiv_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descL = *L;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgetrs_incpiv_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_sgetrs_incpiv_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/
    plasma_parallel_call_6(plasma_pstrsmpl,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_desc, descL,
        const int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pstrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaUpper,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        float, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
