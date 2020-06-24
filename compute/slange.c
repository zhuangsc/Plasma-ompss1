/**
 *
 * @file slange.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:45:08 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_slange returns the value
 *
 *     slange = ( max(abs(A(i,j))), NORM = PlasmaMaxNorm
 *              (
 *              ( norm1(A),         NORM = PlasmaOneNorm
 *              (
 *              ( normI(A),         NORM = PlasmaInfNorm
 *              (
 *              ( normF(A),         NORM = PlasmaFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = PlasmaMaxNorm: Max norm
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *          = PlasmaFrobeniusNorm: Frobenius norm
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0. When M = 0,
 *          the returned value is set to zero.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0. When N = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval the norm described above.
 *
 *******************************************************************************
 *
 * @sa PLASMA_slange_Tile
 * @sa PLASMA_slange_Tile_Async
 * @sa PLASMA_clange
 * @sa PLASMA_dlange
 * @sa PLASMA_slange
 *
 ******************************************************************************/
float PLASMA_slange(PLASMA_enum norm, int M, int N,
                     float *A, int LDA)
{
    int NB;
    int status;
    float value;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slange", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (norm != PlasmaMaxNorm) && (norm != PlasmaOneNorm)
        && (norm != PlasmaInfNorm) && (norm != PlasmaFrobeniusNorm) ) {
        plasma_error("PLASMA_slange", "illegal value of norm");
        return -1;
    }
    if (M < 0) {
        plasma_error("PLASMA_slange", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_slange", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_slange", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (min(N, M) == 0)
      return (float)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_slange", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_siplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_slange_Tile_Async(norm, &descA, &value, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    plasma_sequence_destroy(plasma, sequence);
    return value;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_slange_Tile - Tile equivalent of PLASMA_slange().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = PlasmaMaxNorm: Max norm
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *          = PlasmaFrobeniusNorm: Frobenius norm
 *
 * @param[in] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_slange
 * @sa PLASMA_slange_Tile_Async
 * @sa PLASMA_clange_Tile
 * @sa PLASMA_dlange_Tile
 * @sa PLASMA_slange_Tile
 *
 ******************************************************************************/
float PLASMA_slange_Tile(PLASMA_enum norm, PLASMA_desc *A)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    float value;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slange_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_slange_Tile_Async(norm, A, &value, sequence, &request);
    plasma_dynamic_sync();
    plasma_sequence_destroy(plasma, sequence);
    return value;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_slange_Tile_Async - Non-blocking equivalent of PLASMA_slange_Tile().
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
 * @sa PLASMA_slange
 * @sa PLASMA_slange_Tile
 * @sa PLASMA_clange_Tile_Async
 * @sa PLASMA_dlange_Tile_Async
 * @sa PLASMA_slange_Tile_Async
 *
 ******************************************************************************/
int PLASMA_slange_Tile_Async(PLASMA_enum norm, PLASMA_desc *A, float *value,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    float *work = NULL;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slange_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_slange_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_slange_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_slange_Tile", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_slange_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (norm != PlasmaMaxNorm) && (norm != PlasmaOneNorm)
         && (norm != PlasmaInfNorm) && (norm != PlasmaFrobeniusNorm) ) {
        plasma_error("PLASMA_slange", "illegal value of norm");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (min(descA.m, descA.n) == 0) {
        *value = 0.0;
        return PLASMA_SUCCESS;
    }

    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) {
        if (norm == PlasmaFrobeniusNorm) {
            work = plasma_shared_alloc(plasma, 2*PLASMA_SIZE, PlasmaRealDouble );
        } else {
            work = plasma_shared_alloc(plasma,   PLASMA_SIZE, PlasmaRealDouble );
        }
    }

    plasma_parallel_call_6(plasma_pslange,
        PLASMA_enum, norm,
        PLASMA_desc, descA,
        float*, work,
        float*, value,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    if (work != NULL)
        plasma_shared_free( plasma, work );

    return PLASMA_SUCCESS;
}
