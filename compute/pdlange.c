/**
 *
 * @file pdlange.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:45:11 2014
 *
 **/
#include <stdlib.h>
#include <math.h>
#include "common.h"

#define A(m, n, i, j, ldt)  (BLKADDR(A, double, m, n)+((j)*(ldt)+(i)))
#define Ad(m, n) BLKADDR(A, double, m, n)
#define descx(m, n) BLKADDR(descx, double, m, n)
#define descSx(m, n) BLKADDR(descSx, double, m, n)
/***************************************************************************//**
 *
 **/
void plasma_pdlange(plasma_context_t *plasma)
{
    PLASMA_enum norm;
    PLASMA_desc A;
    double *work;
    double *result;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n;
    int next_m;
    int next_n;
    int ldam;
    int step, lrank;
    int X, X1, X2, Y, Y1, Y2;

    double* lwork;
    double normtmp, normtmp2;
    double *scale, *sumsq;

    plasma_unpack_args_6(norm, A, work, result, sequence, request);
    *result = 0.0;

    if (PLASMA_RANK == 0) {
        if ( norm == PlasmaFrobeniusNorm ) {
            memset(work, 0, 2*PLASMA_SIZE*sizeof(double));
        } else {
            memset(work, 0,   PLASMA_SIZE*sizeof(double));
        }
    }
    ss_init(PLASMA_SIZE, 1, 0);

    switch (norm) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        n = 0;
        m = PLASMA_RANK;
        while (m >= A.mt && n < A.nt) {
            n++;
            m = m-A.mt;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m-A.mt;
            }

            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;

            Y1 = n == 0      ?  A.j       %A.nb   : 0;
            Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y = Y2 - Y1;

            ldam = BLKLDD(A, m);
            CORE_dlange(PlasmaMaxNorm, X, Y, A(m, n, X1, Y1, ldam), ldam, NULL, &normtmp);

            if (normtmp > work[PLASMA_RANK])
                work[PLASMA_RANK] = normtmp;

            m = next_m;
            n = next_n;
        }
        ss_cond_set(PLASMA_RANK, 0, 1);
        break;
    /*
     *  PlasmaOneNorm
     */
    case PlasmaOneNorm:

        n = PLASMA_RANK;
        normtmp2 = 0.0;
        lwork = (double*)plasma_private_alloc(plasma, A.nb, PlasmaRealDouble);

        while (n < A.nt) {
            Y1 = n == 0      ?  A.j       %A.nb   : 0;
            Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y = Y2 - Y1;
            memset(lwork, 0, A.nb*sizeof(double));
            for (m = 0; m < A.mt; m++) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                ldam = BLKLDD(A, m);
                CORE_dasum(
                    PlasmaColumnwise, PlasmaUpperLower,
                    X, Y,
                    A(m, n, X1, Y1, ldam), ldam,
                    lwork);
            }
            CORE_dlange(PlasmaMaxNorm, Y, 1, lwork, 1, NULL, &normtmp);
            if (normtmp > normtmp2)
                normtmp2 = normtmp;

            n += PLASMA_SIZE;
        }
        work[PLASMA_RANK] = normtmp2;
        ss_cond_set(PLASMA_RANK, 0, 1);
        plasma_private_free(plasma, lwork);
        break;
    /*
     *  PlasmaInfNorm
     */
    case PlasmaInfNorm:
        m = PLASMA_RANK;
        normtmp2 = 0.0;
        lwork = (double*)plasma_private_alloc(plasma, A.mb, PlasmaRealDouble);

        while (m < A.mt) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;

            ldam = BLKLDD(A, m);
            memset(lwork, 0, A.mb*sizeof(double));

            for (n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                CORE_dasum(
                    PlasmaRowwise, PlasmaUpperLower,
                    X, Y,
                    A(m, n, X1, Y1, ldam), ldam,
                    lwork);
            }
            CORE_dlange(PlasmaMaxNorm, X, 1, lwork, 1, NULL, &normtmp);
            if (normtmp > normtmp2)
                normtmp2 = normtmp;

            m += PLASMA_SIZE;
        }
        work[PLASMA_RANK] = normtmp2;
        ss_cond_set(PLASMA_RANK, 0, 1);
        plasma_private_free(plasma, lwork);
        break;
    /*
     *  PlasmaFrobeniusNorm
     */
    case PlasmaFrobeniusNorm:
        n = 0;
        m = PLASMA_RANK;
        
        scale = work + 2 * PLASMA_RANK;
        sumsq = work + 2 * PLASMA_RANK + 1;

        *scale = 0.;
        *sumsq = 1.;

        while (m >= A.mt && n < A.nt) {
            n++;
            m = m-A.mt;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m-A.mt;
            }

            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;

            Y1 = n == 0      ?  A.j       %A.nb   : 0;
            Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y = Y2 - Y1;

            ldam = BLKLDD(A, m);
            CORE_dgessq( X, Y, A(m, n, X1, Y1, ldam), ldam, scale, sumsq );

            m = next_m;
            n = next_n;
        }
        ss_cond_set(PLASMA_RANK, 0, 1);
        break;
    default:;
    }

    if (norm != PlasmaFrobeniusNorm) {
        step = 1;
        lrank = PLASMA_RANK;
        while ( (lrank%2 == 0) && (PLASMA_RANK+step < PLASMA_SIZE) ) {
            ss_cond_wait(PLASMA_RANK+step, 0, step);
            work[PLASMA_RANK] = max(work[PLASMA_RANK], work[PLASMA_RANK+step]);
            lrank = lrank >> 1;
            step  = step << 1;
            ss_cond_set(PLASMA_RANK, 0, step);
        }
        if (PLASMA_RANK > 0) {
            while( lrank != 0 ) {
                if (lrank%2 == 1) {
                    ss_cond_set(PLASMA_RANK, 0, step);
                    lrank = 0;
                } else {
                    lrank = lrank >> 1;
                    step  = step << 1;
                    ss_cond_set(PLASMA_RANK, 0, step);
                }
            }
        }

        if (PLASMA_RANK == 0)
            *result = work[0];
    } 
    else {
        step = 1;
        lrank = PLASMA_RANK;
        while ( (lrank%2 == 0) && (PLASMA_RANK+step < PLASMA_SIZE) ) {
            double scale1, scale2;
            double sumsq1, sumsq2;

            ss_cond_wait(PLASMA_RANK+step, 0, step);

            scale1 = work[ 2 * PLASMA_RANK ];
            sumsq1 = work[ 2 * PLASMA_RANK + 1 ];
            scale2 = work[ 2 * (PLASMA_RANK+step) ];
            sumsq2 = work[ 2 * (PLASMA_RANK+step) + 1 ];

            if ( scale2 != 0. ){
                if( scale1 < scale2 ) {
                    work[2 * PLASMA_RANK+1] = sumsq2 + (sumsq1 * (( scale1 / scale2 ) * ( scale1 / scale2 )));
                    work[2 * PLASMA_RANK  ] = scale2;
                } else {
                    work[2 * PLASMA_RANK+1] = sumsq1 + (sumsq2 * (( scale2 / scale1 ) * ( scale2 / scale1 )));
                }
            }
            lrank = lrank >> 1;
            step  = step << 1;
            ss_cond_set(PLASMA_RANK, 0, step);
        }
        if (PLASMA_RANK > 0) {
            while( lrank != 0 ) {
                if (lrank%2 == 1) {
                    ss_cond_set(PLASMA_RANK, 0, step);
                    lrank = 0;
                } else {
                    lrank = lrank >> 1;
                    step  = step << 1;
                    ss_cond_set(PLASMA_RANK, 0, step);
                }
            }
        }

        if (PLASMA_RANK == 0)
            *result = work[0] * sqrt( work[1] );
    }
    ss_finalize();
}

/***************************************************************************//**
 *
 **/
void plasma_pdlange_quark(PLASMA_enum norm, PLASMA_desc A, double *work, double *result,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    double* lwork;
    int X, X1, X2, Y, Y1, Y2;
    int ldam;
    int m, n, k;
    int szeW;
    int nbworker = 1;
    //int PLASMA_SIZE2 = omp_get_num_threads();
    double e0 = 0; double tol = 0.30; int nnz; int IONE=1; int ISEED[4] = {0,0,0,1}; int i, j; double normx; double normSx; int cnt = 0;int maxitr = 100;
    double *res; double beta; double zbeta; 
    double zone = (double)1.0;
    int ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    *result = 0.0;
    switch ( norm ) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        szeW = A.mt*A.nt;
        lwork = (double*)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);
        #pragma omp register ([szeW]lwork)
        memset(lwork, 0, szeW*sizeof(double));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                RT_CORE_dlange_f1(
                    plasma->quark, &task_flags,
                    PlasmaMaxNorm, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    0, &(lwork[A.mt*n+m]),
                    lwork, szeW);
            }
        }
        RT_CORE_dlange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.mt, A.nt,
            lwork, A.mt, szeW,
            0, result);

        RT_CORE_free(plasma->quark, &task_flags, lwork, szeW*sizeof(double));
        break;

    /*
     *  PlasmaOneNorm
     */
    case PlasmaOneNorm:
        lwork = (double*)plasma_shared_alloc(plasma, (A.n+1), PlasmaRealDouble);
        #pragma omp register ([A.n+1]lwork)
        memset(lwork, 0, (A.n+1)*sizeof(double));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                RT_CORE_dasum_f1(
                    plasma->quark, &task_flags,
                    PlasmaColumnwise, PlasmaUpperLower, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    &(lwork[n*A.nb+1]), A.nb,
                    lwork, A.n);
            }
        }
        RT_CORE_dlange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.n+1, 1,
            lwork, 1, A.n+1,
            0, result);

        RT_CORE_free(plasma->quark, &task_flags, lwork, (A.n+1)*sizeof(double));
        break;
    /*
     *  PlasmaSecondNorm
     */
    case PlasmaSecondNorm:
        lwork = (double*)plasma_shared_alloc(plasma, (A.n+1), PlasmaRealDouble);
        #pragma omp register ([A.n+1]lwork)
        memset(lwork, 0, (A.n+1)*sizeof(double));
        double *Sx = (double*)plasma_shared_alloc(plasma, (A.n+1), PlasmaRealDouble);
        #pragma omp register ([A.n+1]Sx)
        PLASMA_desc *descx2 ;  PLASMA_Desc_Create( &descx2 , lwork, PlasmaRealDouble, A.nb, 1, A.nb, A.n, 1, 0, 0, A.n, 1);
        PLASMA_desc descx = plasma_desc_submatrix(*descx2, 0, 0, A.n, 1);
        PLASMA_desc *descSx2;  PLASMA_Desc_Create( &descSx2, Sx   , PlasmaRealDouble, A.nb, 1, A.nb, A.n, 1, 0, 0, A.n, 1);
        PLASMA_desc descSx = plasma_desc_submatrix(*descSx2, 0, 0, A.n, 1);
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                RT_CORE_dasum_f1(
                    plasma->quark, &task_flags,
                    PlasmaColumnwise, PlasmaUpperLower, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    &(lwork[n*A.nb+1]), A.nb,
                    lwork, A.n);
            }          
        }
        RT_CORE_dlange(
            plasma->quark, &task_flags,
            PlasmaFrobeniusNorm, A.n+1, 1,
            lwork, 1, A.n+1,
            0, result); // Why in One and Inf norm take the whole lwork
        RT_CORE_dscal( plasma->quark, &task_flags, A.n, result, lwork, 1);
        beta = 0.0; double conv = result[0]; double tolconv = tol*result[0];
        RT_CORE_tolconv(plasma->quark, &task_flags, &tol, result, &tolconv, &conv);
        //while ( conv > tolconv ){ 
        while ( cnt < 2 )
        //#pragma omp task in(&conv)
        {
        //if(conv < tolconv){return;}
         
           RT_CORE_dgemm_2norm( plasma->quark, &task_flags, PlasmaNoTrans, PlasmaNoTrans, 1.0, A, lwork, A.n, 0.0, Sx, A.n);
           RT_CORE_dlange(
               plasma->quark, &task_flags,
               PlasmaFrobeniusNorm, A.n, 1,
               lwork, 1, A.n,
               0, &normx);
           RT_CORE_dlange(
               plasma->quark, &task_flags,
               PlasmaFrobeniusNorm, A.n, 1,
               Sx, 1, A.n,
               0, &normSx);
           RT_CORE_dscal(plasma->quark, &task_flags, A.n, &normx, lwork, 1);
           //RT_CORE_dconv( &normx, &normSx, result, &e0);
           RT_CORE_dconv(plasma->quark, &task_flags, &normx, &normSx, result, &e0);
            
           RT_CORE_conv(plasma->quark, &task_flags, result, &e0, &conv);
           cnt = cnt + 1;
           if ( cnt > maxitr){
               printf("\n normest:notconverge \n");
               return;
           }
           //#pragma omp atomic
           //RT_CORE_conv(plasma->quark, &task_flags, result, &e0, &conv);
       }

       RT_CORE_free(plasma->quark, &task_flags, lwork, (A.n+1)*sizeof(double));
       RT_CORE_free(plasma->quark, &task_flags, Sx, (A.n+1)*sizeof(double));
RT_dynamic_sync();
       break;
    /*
     *  PlasmaInfNorm
     */
    case PlasmaInfNorm:
        lwork = (double*)plasma_shared_alloc(plasma, (A.m+1), PlasmaRealDouble);
        #pragma omp register ([A.m+1]lwork)
        memset(lwork, 0, (A.m+1)*sizeof(double));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                RT_CORE_dasum_f1(
                    plasma->quark, &task_flags,
                    PlasmaRowwise, PlasmaUpperLower, X, Y,
                    A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                    &(lwork[m*A.mb+1]), A.mb,
                    lwork, A.m);
            }
        }
        RT_CORE_dlange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.m+1, 1,
            lwork, 1, A.m+1,
            0, result);

        RT_CORE_free(plasma->quark, &task_flags, lwork, (A.m+1)*sizeof(double));
        break;
    /*
     *  PlasmaFrobeniusNorm
     */
    case PlasmaFrobeniusNorm:
        szeW = 2*(PLASMA_SIZE+1);
        printf("\n PLASMA_SIZE %d \n", PLASMA_SIZE);
        lwork = (double*)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);
        #pragma omp register ([szeW]lwork)

        for(m = 0; m < PLASMA_SIZE+1; m++) {
            lwork[2*m  ] = 0.;
            lwork[2*m+1] = 1.;
        }

        k = 0;
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            for(n = 0; n < A.nt; n++) {
                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;
                k++; nbworker++;
                RT_CORE_dgessq_f1(
                    plasma->quark, &task_flags,
                    X, Y,
                    A(m, n, X1, Y1, ldam), ldam,
                    lwork + 2*k, 
                    lwork + 2*k + 1,
                    lwork, szeW, OUTPUT | GATHERV );
                k = k % PLASMA_SIZE;
            }
        }
        RT_CORE_dplssq(
            plasma->quark, &task_flags,
            min(nbworker, PLASMA_SIZE+1), lwork, result );

        RT_CORE_free(plasma->quark, &task_flags, lwork, szeW*sizeof(double));
        break;
    default:;
    }
}
