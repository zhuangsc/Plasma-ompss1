/**
 *
 * @file psbarrier.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * Barrier for algorithm mixing computation on tile/panel.
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2009-11-15
 *
 * @generated s Tue Jan  7 11:45:13 2014
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)

/***************************************************************************//**
 *  Barrier from tiles to panels
 **/
void plasma_psbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int m, n;
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
 
    for (n = 0; n < A.nt; n++)
    {
        /* Protection from previous GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(0, n), INOUT,
                          0);

        for (m = 1; m < A.mt; m++)
        {
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(float)*A.mb*A.nb, A(0, n), INOUT | GATHERV,
                              sizeof(float)*A.mb*A.nb, A(m, n), INOUT,
                              0);
        }

        /* Protection to next GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(0, n), INOUT,
                          0);
    }
}

/***************************************************************************//**
 *  Barrier from panels to tiles
 **/
void plasma_psbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int m, n;
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
 
    for (n = 0; n < A.nt; n++)
    {
        /* Protection from previous GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(0, n), INOUT,
                          0);

        for (m = 1; m < A.mt; m++)
        {
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(float)*A.mb*A.nb, A(0, n), INPUT,
                              sizeof(float)*A.mb*A.nb, A(m, n), INOUT,
                              0);
        }
    }
}

/***************************************************************************//**
 *  Barrier from tiles to panels
 **/
void plasma_psbarrier_tl2row_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int m, n;
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
 
    for (m = 0; m < A.mt; m++)
    {
        /* Protection from previous GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(m, 0), INOUT,
                          0);

        for (n = 1; n < A.nt; n++)
        {
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(float)*A.mb*A.nb, A(m, 0), INOUT | GATHERV,
                              sizeof(float)*A.mb*A.nb, A(m, n), INOUT,
                              0);
        }

        /* Protection to next GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(m, 0), INOUT,
                          0);
    }
}

/***************************************************************************//**
 *  Barrier from panels to tiles
 **/
void plasma_psbarrier_row2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int m, n;
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
 
    for (m = 0; m < A.mt; m++)
    {
        /* Protection from previous GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(m, 0), INOUT,
                          0);

        for (n = 1; n < A.nt; n++)
        {
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(float)*A.mb*A.nb, A(m, 0), INPUT,
                              sizeof(float)*A.mb*A.nb, A(m, n), INOUT,
                              0);
        }
    }
}

