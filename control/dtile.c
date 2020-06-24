/**
 *
 * @file dtile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:45:15 2014
 *
 **/
#include "common.h"
#include "auxiliary.h"
#include "tile.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *          If PLASMA_TRANSLATION_MODE is set to PLASMA_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          PLASMA_TRANSLATION_MODE is set to PLASMA_OUTOFPLACE,
 *          A->mat has to be allocated before.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dLapack_to_Tile_Async
 * @sa PLASMA_dTile_to_Lapack
 * @sa PLASMA_cLapack_to_Tile
 * @sa PLASMA_dLapack_to_Tile
 * @sa PLASMA_sLapack_to_Tile
 *
 ******************************************************************************/
int PLASMA_dLapack_to_Tile(double *Af77, int LDA, PLASMA_desc *A)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dLapack_to_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dLapack_to_Tile", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    plasma_sequence_create(plasma, &sequence);
    plasma_parallel_call_5(
        plasma_pdlapack_to_tile,
        double*, Af77,
        int, LDA,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dLapack_to_Tile_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of PLASMA_dLapack_to_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *          If PLASMA_TRANSLATION_MODE is set to PLASMA_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          PLASMA_TRANSLATION_MODE is set to PLASMA_OUTOFPLACE,
 *          A->mat has to be allocated before.
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
 * @sa PLASMA_dTile_to_Lapack_Async
 * @sa PLASMA_dLapack_to_Tile
 * @sa PLASMA_cLapack_to_Tile_Async
 * @sa PLASMA_dLapack_to_Tile_Async
 * @sa PLASMA_sLapack_to_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dLapack_to_Tile_Async(double *Af77, int LDA, PLASMA_desc *A,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dLapack_to_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dLapack_to_Tile", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    plasma_parallel_call_5(
        plasma_pdlapack_to_tile,
        double*, Af77,
        int, LDA,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If PLASMA_TRANSLATION_MODE is set to PLASMA_INPLACE,
 *          Af77 has to be A->mat, else if
 *          PLASMA_TRANSLATION_MODE is set to PLASMA_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dTile_to_Lapack_Async
 * @sa PLASMA_dLapack_to_Tile
 * @sa PLASMA_cTile_to_Lapack
 * @sa PLASMA_dTile_to_Lapack
 * @sa PLASMA_sTile_to_Lapack
 *
******************************************************************************/
int PLASMA_dTile_to_Lapack(PLASMA_desc *A, double *Af77, int LDA)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dTile_to_Lapack", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dTile_to_Lapack", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    plasma_sequence_create(plasma, &sequence);
    plasma_static_call_5(
        plasma_pdtile_to_lapack,
        PLASMA_desc, descA,
        double*, Af77,
        int, LDA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dTile_to_Lapack_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of PLASMA_dTile_to_Lapack().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If PLASMA_TRANSLATION_MODE is set to PLASMA_INPLACE,
 *          Af77 has to be A->mat, else if
 *          PLASMA_TRANSLATION_MODE is set to PLASMA_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
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
 * @sa PLASMA_dLapack_to_Tile_Async
 * @sa PLASMA_dTile_to_Lapack
 * @sa PLASMA_cTile_to_Lapack_Async
 * @sa PLASMA_dTile_to_Lapack_Async
 * @sa PLASMA_sTile_to_Lapack_Async
 *
 ******************************************************************************/
int PLASMA_dTile_to_Lapack_Async(PLASMA_desc *A, double *Af77, int LDA,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dTile_to_Lapack", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dTile_to_Lapack", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    plasma_static_call_5(
        plasma_pdtile_to_lapack,
        PLASMA_desc, descA,
        double*, Af77,
        int, LDA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
