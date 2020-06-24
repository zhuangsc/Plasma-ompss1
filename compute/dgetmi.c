/**
 *
 * @file dgetmi.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated d Tue Jan  7 11:45:09 2014
 *
 **/

#include <stdlib.h>
#include <sys/types.h>
#include "common.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  PLASMA_dgetmi Implementation of inplace transposition
 *    based on the GKK algorithm by Gustavson, Karlsson, Kagstrom.
 *    This algorithm shift some cycles to transpose the matrix.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size L*m*n.
 *
 * @param[in] f_in
 *         Original format of the matrix A. Must be part of (PlasmaCM, PlasmaRM,
 *         PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB)
 *
 * @param[in] imb
 *         Number of rows of the problem
 *
 * @param[in] inb
 *         Number of columns in the problem
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgetmi_Async
 *
 ******************************************************************************/
int PLASMA_dgetmi(int m, int n, double *A, PLASMA_enum f_in, int imb, int inb)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error(__func__, "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    plasma_sequence_create(plasma, &sequence);

    PLASMA_dgetmi_Async( m, n, A, 
                         f_in, imb, inb,
                         sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  PLASMA_dgetmi_Async Implementation of inplace transposition
 *    based on the GKK algorithm by Gustavson, Karlsson, Kagstrom.
 *    This algorithm shift some cycles to transpose the matrix.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size L*m*n.
 *
 * @param[in] f_in
 *         Original format of the matrix A. Must be part of (PlasmaCM, PlasmaRM,
 *         PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB)
 *
 * @param[in] mb
 *         Number of rows of the problem
 *
 * @param[in] nb
 *         Number of columns in the problem
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
 * @sa PLASMA_dgetmi
 *
 ******************************************************************************/
int PLASMA_dgetmi_Async(int m, int n, double *A, PLASMA_enum f_in, int mb, int nb, 
                        PLASMA_sequence *sequence, PLASMA_request *request)
{
    /* convert */
    switch ( f_in ) {
    case PlasmaCM   :
        PLASMA_dgecfi_Async(m, n, A, PlasmaCM,   mb, nb, PlasmaRM,   nb, mb, sequence, request);
        break;
    case PlasmaCCRB :
        PLASMA_dgecfi_Async(m, n, A, PlasmaCCRB, mb, nb, PlasmaRRRB, nb, mb, sequence, request);
        break;
    case PlasmaCRRB :
        PLASMA_dgecfi_Async(m, n, A, PlasmaCRRB, mb, nb, PlasmaRCRB, nb, mb, sequence, request);
        break;
    case PlasmaRCRB :
        PLASMA_dgecfi_Async(m, n, A, PlasmaRCRB, mb, nb, PlasmaCRRB, nb, mb, sequence, request);
        break;
    case PlasmaRRRB :
        PLASMA_dgecfi_Async(m, n, A, PlasmaRRRB, mb, nb, PlasmaCCRB, nb, mb, sequence, request);
        break;
    case PlasmaRM   :
        PLASMA_dgecfi_Async(m, n, A, PlasmaRM,   mb, nb, PlasmaCM,   nb, mb, sequence, request);
        break;
    default:
        plasma_error(__func__, "unknown format");
    }
    return PLASMA_SUCCESS;
}
