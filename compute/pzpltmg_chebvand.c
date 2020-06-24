/**
 *
 * @file pzpltmg_chebvand.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Ichitaro Yamazaki
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include <stdlib.h>
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)

/***************************************************************************//**
 *  Parallel tile chebvandion matrix generation -- Dynamic scheduling
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999859
 *
 *  Vandermonde-like matrix for the Chebyshev polynomials
 *
 *  Produces the (primal) Chebyshev Vandermonde matrix based on the vector of
 *  points p, which define where the Chebyshev polynomial is calculated.
 *
 *  If seed != 0, C(i,j) = Ti – 1(p(j)) where Ti – 1 is the Chebyshev
 *  polynomial of degree i – 1, and p is a vector of N equally spaced points on
 *  the interval [0,1].
 *
 */
void plasma_pzpltmg_chebvand_quark( PLASMA_desc A,
                                   PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex64_t **W;
    int m, n;
    int ldam;
    int tempm0, tempn0;
    int tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    W = (PLASMA_Complex64_t**) malloc ( A.nt * sizeof( PLASMA_Complex64_t* ) );

    /* Initialize the full matrix */
    for (n = 0; n < A.nt; n++) {
        tempn0 = n * A.nb;
        tempnn = n == A.nt-1 ? A.n - tempn0 : A.nb;

        W[n] = (PLASMA_Complex64_t*)plasma_shared_alloc(plasma, 2*tempnn, PlasmaComplexDouble);

        for (m = 0; m < A.mt; m++) {
            tempm0 = m * A.mb;
            tempmm = m == A.mt-1 ? A.m - tempm0 : A.mb;
            ldam = BLKLDD(A, m);

            QUARK_CORE_zpltmg_chebvand(
                plasma->quark, &task_flags,
                tempmm, tempnn, A(m, n), ldam,
                A.n, tempm0, tempn0, W[n] );
        }

        QUARK_CORE_free(plasma->quark, &task_flags,
                        W[n], (2*tempnn)*sizeof(PLASMA_Complex64_t));
    }

    free(W);
}
