/**
 *
 * @file dgetrf_nopiv.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @date 2013-02-01
 *
 * @generated d Tue Jan  7 11:45:10 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dgetrf_nopiv - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm without row pivoting.
 *  WARNING: Don't use this function if you are not sure your matrix is diagonal
 *  dominant.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been
 *               completed, but the factor U is exactly singular, and division
 *               by zero will occur if it is used to solve a system of
 *               equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgetrf_nopiv_Tile
 * @sa PLASMA_dgetrf_nopiv_Tile_Async
 * @sa PLASMA_cgetrf_nopiv
 * @sa PLASMA_dgetrf_nopiv
 * @sa PLASMA_sgetrf_nopiv
 *
 ******************************************************************************/
int PLASMA_dgetrf_nopiv(int M, int N,
                        double *A, int LDA)
{
    int NB;
    int status;
    PLASMA_desc descA ;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_nopiv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_dgetrf_nopiv", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_dgetrf_nopiv", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_dgetrf_nopiv", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGESV, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgetrf_nopiv", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                             sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_diplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, N, N,
                             sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_dgetrf_nopiv_Tile_Async(&descA, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,
                             sequence, &request );
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N,
                             sequence, &request );
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dgetrf_nopiv_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of PLASMA_dgetrf_nopiv().  Operates on matrices stored by
 *  tiles.  All matrices are passed through descriptors.  All dimensions are
 *  taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgetrf_nopiv
 * @sa PLASMA_dgetrf_nopiv_Tile_Async
 * @sa PLASMA_cgetrf_nopiv_Tile
 * @sa PLASMA_dgetrf_nopiv_Tile
 * @sa PLASMA_sgetrf_nopiv_Tile
 * @sa PLASMA_dgetrs_Tile
 *
 ******************************************************************************/
int PLASMA_dgetrf_nopiv_Tile(PLASMA_desc *A)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_nopiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dgetrf_nopiv_Tile_Async(A, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dgetrf_nopiv_Tile_Async - Computes the tile LU factorization of a
 *  matrix.  Non-blocking equivalent of PLASMA_dgetrf_nopiv_Tile().  May return
 *  before the computation is finished.  Allows for pipelining of operations ar
 *  runtime.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
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
 * @sa PLASMA_dgetrf_nopiv
 * @sa PLASMA_dgetrf_nopiv_Tile
 * @sa PLASMA_cgetrf_nopiv_Tile_Async
 * @sa PLASMA_dgetrf_nopiv_Tile_Async
 * @sa PLASMA_sgetrf_nopiv_Tile_Async
 * @sa PLASMA_dgetrs_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dgetrf_nopiv_Tile_Async(PLASMA_desc *A,
                                   PLASMA_sequence *sequence,
                                   PLASMA_request *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_nopiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_nopiv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dgetrf_nopiv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgetrf_nopiv_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }

    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dgetrf_nopiv_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    plasma_dynamic_call_3(plasma_pdgetrf_nopiv,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
