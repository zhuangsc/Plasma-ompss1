/**
 *
 * @file pzgbrdb.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <lapacke.h>

#undef REAL
#define COMPLEX

#define DEP(m)  &(DEP[m])
#define A(_m, _n) (PLASMA_Complex64_t *)plasma_geteltaddr(&A, (_m), (_n), eltsize)
/***************************************************************************//**
 *  Parallel Reduction from BAND Bidiagonal to the final condensed form - dynamic scheduler
 **/
void plasma_pzgbrdb_quark(PLASMA_enum uplo,
                          PLASMA_desc A, double *D, double *E, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

#ifdef COMPLEX
    static double dzero = (double) 0.0;
    double absztmp;
#endif

    static PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    static PLASMA_Complex64_t zzero = (PLASMA_Complex64_t) 0.0;
    PLASMA_Complex64_t *VQ, *TAUQ,*VP, *TAUP;
    PLASMA_Complex64_t ztmp, V, TAU;
    int M, N, NB, MINMN, INgrsiz, INthgrsiz, BAND;
    int myid, grsiz, shift=3, stt, st, ed, stind, edind;
    int blklastind, colpt, PCOL, ACOL, MCOL;
    int stepercol, mylastid, grnb, grid;
    int *DEP,*MAXID;
    int i, sweepid, m;
    int thgrsiz, thgrnb, thgrid, thed;
    size_t eltsize = plasma_element_size(A.dtyp);

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    M     = A.m;
    N     = A.n;
    NB    = A.mb;
    MINMN = min(M,N);

    /* Quick return */
    if ( MINMN == 0 ){
        return;
    }
    if ( NB == 0 ) {
        memset(D, 0,  MINMN   *sizeof(double));
        memset(E, 0, (MINMN-1)*sizeof(double));
#ifdef COMPLEX
        for (i=0; i<MINMN; i++)
            D[i] = cabs(*A(i,i));
#else
        for (i=0; i<MINMN; i++)
            D[i] = *A(i,i);
#endif
        return;
    }

    /*
     * Barrier is used because the bulge have to wait until
     * the reduction to band has been finish.
     * otherwise, I can remove this BARRIER when I integrate
     * the function dependencies link inside the reduction to
     * band. Keep in mind the case when NB=1, where no bulge-chasing.
     */
    /***************************************************************/
    QUARK_Barrier(plasma->quark);
    /***************************************************************/

    /*
     * Case NB=1 ==> matrix is already Bidiagonal. no need to bulge.
     * Make diagonal and superdiagonal elements real, storing them in
     * D and E. if PlasmaLower, first transform lower bidiagonal form
     * to upper bidiagonal by applying plane rotations/ Householder
     * from the left, overwriting superdiagonal elements then make
     * elements real of the resulting upper Bidiagonal. if PlasmaUpper
     * then make its elements real.  For Q, PT: ZSCAL should be done
     * in case of WANTQ.
     */
    if ( NB == 1 ) {
        memset(D, 0, MINMN*sizeof(double));
        memset(E, 0, (MINMN-1)*sizeof(double));

        if(uplo==PlasmaLower){
            for (i=0; i<(MINMN-1); i++)
            {
                /* generate Householder to annihilate a(i+1,i) and create a(i,i+1) */
                V             = *A((i+1), i);
                *A((i+1),  i) = zzero;
                LAPACKE_zlarfg_work( 2, A(i, i), &V, 1, &TAU);
                /* apply Left*/
                TAU  = conj(TAU);
                ztmp = TAU*V;
                V    = conj(V);
                *A(i,   i+1) = - V * TAU * (*A(i+1, i+1));
                *A(i+1, i+1) = *(A(i+1, i+1)) * (zone - V * ztmp);
            }
        }
        /* PlasmaLower or PlasmaUpper, both are now upper */
        /* Make diagonal and superdiagonal elements real,
         * storing them in D and E
         */
#ifdef COMPLEX
        ztmp = zone;
        for (i=0; i<MINMN; i++)
        {
            ztmp     = *A(i, i) * conj(ztmp);
            absztmp  = cabs(ztmp);
            D[i]     = absztmp;               /* diag value */
            if(absztmp != dzero)
                ztmp = (PLASMA_Complex64_t) (ztmp / absztmp);
            else
                ztmp = zone;
            if(i<(MINMN-1)) {
                ztmp     = *A(i, (i+1)) * conj(ztmp);
                absztmp  = cabs(ztmp);
                E[i]     = absztmp;            /* upper off-diag value */
                if(absztmp != dzero)
                    ztmp = (PLASMA_Complex64_t) (ztmp / absztmp);
                else
                    ztmp = zone;
            }
        }
#else
        for (i=0; i < MINMN-1; i++) {
            D[i] = *A(i, i  );
            E[i] = *A(i, i+1);
        }
        D[i] = *A(i, i);
#endif
        return;
    }

    /*
     * Case MINMN<NB ==> matrix is very small and better to call lapack ZGETRD.
     *
     *    Use fact that one row of block is stored the same way than in LAPACK
     *    Doesn't work if M > NB because of tile storage
     */
    if ( MINMN <= 0 )
    {
        PLASMA_Complex64_t *work, *taup, *tauq;
        int info, ldwork = N*N;
        work = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, ldwork, PlasmaComplexDouble);
        taup = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, MINMN,  PlasmaComplexDouble);
        tauq = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, MINMN,  PlasmaComplexDouble);
        
        info = LAPACKE_zgebrd_work(LAPACK_COL_MAJOR, M, N,
                                   A(0,0), A.lm, D, E, taup, tauq, work, ldwork);
        plasma_shared_free(plasma, (void*) work);
        plasma_shared_free(plasma, (void*) taup);
        plasma_shared_free(plasma, (void*) tauq);

        if( info == 0 )
            sequence->status = PLASMA_SUCCESS;
        else
            plasma_sequence_flush(plasma->quark, sequence, request, info);
        return;
    }

    /* General case NB > 1 && N > NB */
    DEP   = (int *)                plasma_shared_alloc(plasma, MINMN+1, PlasmaInteger      );
    MAXID = (int *)                plasma_shared_alloc(plasma, MINMN+1, PlasmaInteger      );
    memset(MAXID,0,(MINMN+1)*sizeof(int));
    VQ    = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*MINMN,   PlasmaComplexDouble);
    TAUQ  = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*MINMN,   PlasmaComplexDouble);
    VP    = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*MINMN,   PlasmaComplexDouble);
    TAUP  = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*MINMN,   PlasmaComplexDouble);
    memset(VQ,0,2*MINMN*sizeof(PLASMA_Complex64_t));
    memset(TAUQ,0,2*MINMN*sizeof(PLASMA_Complex64_t));
    memset(VP,0,2*MINMN*sizeof(PLASMA_Complex64_t));
    memset(TAUP,0,2*MINMN*sizeof(PLASMA_Complex64_t));

    /***************************************************************************
     *                       START BULGE CHASING CODE
     **************************************************************************/
    printf("converting tile to lap A.lm %d   ln %d \n",A.lm,A.ln);
    PLASMA_Complex64_t *AB =   (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, A.lm*A.ln,   PlasmaComplexDouble);
    memset(AB,0,A.lm*A.ln*sizeof(PLASMA_Complex64_t));
    plasma_zooptile2lap( A,   AB, NB, NB,  A.lm, A.ln,  sequence, request);
    int LDAB = A.lm;
/*
    plasma_dynamic_call_6( plasma_pzhbcpy_t2bl,
        PLASMA_enum, uplo,
        PLASMA_desc, A,
        PLASMA_Complex64_t*, &(AB[NB]),
        int, M,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);
*/

    QUARK_Barrier(plasma->quark);

    /* 
     * Initialisation of local parameter. those parameter should be
     * input or tuned parameter.
     */
    INgrsiz = 1;
    if( NB > 160 ) {
        INgrsiz = 2;
    }
    else if( NB > 100 ) {
        if( MINMN < 5000 )
            INgrsiz = 2;
        else
            INgrsiz = 4;
    } else {
        INgrsiz = 6;
    }
    INthgrsiz = MINMN;
    BAND      = 0;

    grsiz   = INgrsiz;
    thgrsiz = INthgrsiz;
    if( grsiz   == 0 ) grsiz   = 6;
    if( thgrsiz == 0 ) thgrsiz = MINMN;

    i = shift/grsiz;
    stepercol =  i*grsiz == shift ? i:i+1;

    i       = (MINMN-1)/thgrsiz;
    thgrnb  = i*thgrsiz == (MINMN-1) ? i:i+1;

    printf("starting code grsiz %d thgr %d\n", grsiz, thgrsiz);
    for (thgrid = 1; thgrid<=thgrnb; thgrid++){
        stt  = (thgrid-1)*thgrsiz+1;
        thed = min( (stt + thgrsiz -1), (MINMN-1));
        for (i = stt; i <= MINMN-1; i++){
            ed=min(i,thed);
            if(stt>ed)break;
            for (m = 1; m <=stepercol; m++){
                st=stt;
                for (sweepid = st; sweepid <=ed; sweepid++){
                    /* PCOL:  dependency on the ID of the master of the group of the previous column.  (Previous Column:PCOL). */
                    /* ACOL:  dependency on the ID of the master of the previous group of my column.   (Acctual  Column:ACOL). (it is 0(NULL) for myid=1) */
                    /* MCOL:  OUTPUT dependency on the my ID, to be used by the next ID. (My Column: MCOL). I am the master of this group. */
                    myid     = (i-sweepid)*(stepercol*grsiz) +(m-1)*grsiz + 1;
                    mylastid = myid+grsiz-1;
                    PCOL     = mylastid+shift-1;  /* to know the dependent ID of the previous column. need to know the master of its group*/
                    MAXID[sweepid] = myid;
                    PCOL     = min(PCOL,MAXID[sweepid-1]); /* for the last columns, we might do only 1 or 2 kernel, so the PCOL will be wrong. this is to force it to the last ID of the previous col.*/
                    grnb     = PCOL/grsiz;
                    grid     = grnb*grsiz == PCOL ? grnb:grnb+1;
                    PCOL     = (grid-1)*grsiz +1; /* give me the ID of the master of the group of the previous column.*/
                    ACOL     = myid-grsiz;
                    if(myid==1)ACOL=0;
                    MCOL     = myid;

                    QUARK_CORE_zbrdalg2(
                        plasma->quark, &task_flags,
                        uplo, MINMN, NB,
                        AB, LDAB, VQ, TAUQ, VP, TAUP, i, sweepid, m, grsiz, BAND,
                        DEP(PCOL), DEP(ACOL), DEP(MCOL) );

                    if(mylastid%2 ==0){
                        blklastind      = (mylastid/2)*NB+1+sweepid-1;
                    }else{
                        colpt      = ((mylastid+1)/2)*NB + 1 +sweepid -1 ;
                        stind      = colpt-NB+1;
                        edind      = min(colpt,MINMN);
                        if( (stind>=edind-1) && (edind==MINMN) )
                            blklastind=MINMN;
                        else
                            blklastind=0;
                    }
                    if(blklastind >= (MINMN-1))  stt=stt+1;
                } /* END for sweepid=st:ed    */
            } /* END for m=1:stepercol */
        } /* END for i=1:MINMN-2      */
    } /* END for thgrid=1:thgrnb     */

     /*
     * Barrier used only for now, to be sure that everything
     * is done before copying the D and E and free workspace.
     * this will be removed later when D and E are directly filled
     * during the bulge process.
     */
    QUARK_Barrier(plasma->quark);
    printf("converting back to tile\n");
    
    plasma_zooplap2tile( A, AB, NB, NB, A.lm, A.ln, 0, 0, M, N, sequence, request, plasma_desc_mat_free(&A));
/*
    A = plasma_desc_init(PlasmaComplexDouble, NB,NB, NB*NB, M, N, 0, 0, M, N);
        if ( plasma_desc_mat_alloc( &(A) ) ) {                         
        plasma_error( __func__, "plasma_shared_alloc() failed");        
        {free;};                                                        
        return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        */
    /*
    plasma_parallel_call_5( plasma_pzlapack_to_tile,   
        PLASMA_Complex64_t*, AB,                                      
        int,                 MINMN,                                      
        PLASMA_desc,         A,                                   
        PLASMA_sequence*,    sequence,                                     
        PLASMA_request*,     request);

*/

    QUARK_Barrier(plasma->quark);
    printf("converting done\n");

    plasma_shared_free(plasma, (void*) DEP);
    plasma_shared_free(plasma, (void*) MAXID);
    plasma_shared_free(plasma, (void*) VQ);
    plasma_shared_free(plasma, (void*) TAUQ);
    plasma_shared_free(plasma, (void*) VP);
    plasma_shared_free(plasma, (void*) TAUP);

    /*
     * STORE THE RESULTING diagonal/off-diagonal in D AND E
     */
    memset(D, 0, MINMN*sizeof(double));
    memset(E, 0, (MINMN-1)*sizeof(double));


    /* Make diagonal and superdiagonal elements real,
     * storing them in D and E
     */
    if(uplo==PlasmaUpper){
        for (i=0; i < MINMN-1; i++) {
            D[i] = *A(i, i  );
            E[i] = *A(i, i+1);
        }
        D[i] = *A(i, i);
    }
    else{
        for (i=0; i < MINMN-1; i++) {
            D[i] = *A(i, i  );
            E[i] = *A(i+1, i);
        }
        D[i] = *A(i, i);
    }


} /* END FUNCTION */
