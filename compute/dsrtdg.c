/**
 *
 * @file dlaset.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dsrtdg_Tile_Async 
 * dsrtdg sort the numbers on the diagonal of the matrix
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
 *******************************************************************************/

int PLASMA_dsrtdg_Tile_Async( PLASMA_desc *A, double *work,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsrtdg_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dsrtdg_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dsrtdg_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsrtdg_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dsrtdg_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (min(descA.m, descA.n) == 0) {
        return PLASMA_SUCCESS;
    }

    plasma_dynamic_call_4(plasma_pdsrtdg,
        PLASMA_desc, descA,
        double *, work,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
