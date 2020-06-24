/**
 * @file zheev.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#if defined(USEMKL)
#include <mkl_service.h>
#endif

#include <lapacke.h>
#include "common.h"


/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zheev - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex Hermitian matrix A. The matrix A is
 *  preliminary reduced to tridiagonal form using a two-stage
 *  approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *  Note: Only PlasmaNoVec supported!
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *          Note: Only PlasmaNoVec supported!
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the lower triangle (if uplo = PlasmaLower) or the
 *          upper triangle (if uplo = PlasmaUpper) of A, including the
 *          diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_zheev
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zheev_Tile
 * @sa PLASMA_zheev_Tile_Async
 * @sa PLASMA_cheev
 * @sa PLASMA_dsyev
 * @sa PLASMA_ssyev
 *
 ******************************************************************************/
int PLASMA_zheev(PLASMA_enum jobz, PLASMA_enum uplo, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 double *W,
                 PLASMA_desc *descT,
                 PLASMA_Complex64_t *Q, int LDQ)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descQ;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_zheev", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_zheev", "illegal value of jobz");
        return -1;
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_zheev", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_zheev", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zheev", "illegal value of LDA");
        return -5;
    }
    if (LDQ < max(1, N)) {
        plasma_error("PLASMA_zheev", "illegal value of LDQ");
        return -9;
    }

    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    /*
    if (jobz == PlasmaVec) {
        plasma_error("PLASMA_zheev", "computing the eigenvectors is not supported in this version");
        return -1;
    }
    */

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZHEEV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zheev", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        /*
        if (jobz == PlasmaVec) {
            // No need for conversion, it's just output
            plasma_zdesc_alloc( descQ, NB, NB, LDQ, N, 0, 0, N, N,
                                plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descQ)) );
        }*/
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
        /*
        if (jobz == PlasmaVec) {
            // No need for conversion, it's just output
            descQ = plasma_desc_init(
                PlasmaComplexDouble, NB, NB, NB*NB,
                LDQ, N, 0, 0, N, N);
            descQ.mat = Q;
        }*/
    }

    /* Call the tile interface */
    PLASMA_zheev_Tile_Async(jobz, uplo, &descA, W, descT, Q, LDQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        /*
        if (jobz == PlasmaVec) {
           plasma_zooptile2lap( descQ, Q, NB, NB, LDQ, N,  sequence, &request);
        }*/
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        /*
        if (jobz == PlasmaVec)
           plasma_desc_mat_free(&descQ);*/
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        /*
        if (jobz == PlasmaVec)
           plasma_ziptile2lap( descQ, Q, NB, NB, LDQ, N,  sequence, &request);*/
        plasma_dynamic_sync();
    }


/*
    int i,j;
    int NN=descA.lm;
    FILE *trace_file;
    trace_file = fopen("AJETE/QQ", "w");
    for (j = 0; j < NN ; j++)
          for (i = 0; i < NN ; i++)
                         fprintf(trace_file,"%10d%10d%40.30e\n",i+1,j+1,A[j*NN+i]);
    fclose(trace_file);
*/



    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zheev_Tile - Computes all eigenvalues and, optionally, eigenvectors of a
 *  complex Hermitian matrix A using a two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *          PlasmaVec is NOT supported.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] T
 *          On exit, auxiliary factorization data.
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          PLASMA_Alloc_Workspace_zheev
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zheev_Tile
 * @sa PLASMA_zheev_Tile_Async
 * @sa PLASMA_cheev
 * @sa PLASMA_dsyev
 * @sa PLASMA_ssyev
 *
 ******************************************************************************/
int PLASMA_zheev_Tile(PLASMA_enum jobz, PLASMA_enum uplo,
                      PLASMA_desc *A, double *W,
                      PLASMA_desc *T, PLASMA_Complex64_t *Q, int LDQ)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zheev_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zheev_Tile_Async(jobz, uplo, A, W, T, Q, LDQ, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zheev_Tile_Async - Computes all eigenvalues and, optionally, eigenvectors of a
 *  complex Hermitian matrix A using a two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
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
 * @sa PLASMA_zheev
 * @sa PLASMA_zheev_Tile
 * @sa PLASMA_cheev_Tile_Async
 * @sa PLASMA_dsyev_Tile_Async
 * @sa PLASMA_ssyev_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zheev_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo,
                            PLASMA_desc *A,
                            double *W,
                            PLASMA_desc *T,
                            PLASMA_Complex64_t *Q, int LDQ,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA = *A;
    PLASMA_desc descT = *T;
    PLASMA_desc descQ;
    PLASMA_Complex64_t *AB;
    double *E;
    int i,j;
    int N       = descA.m;
    int NB      = descA.mb;
    int LDAB    = 2*NB+1;
    int stat;
    plasma = plasma_context_self();
    int THREADS = plasma->world_size;
    //int THREADS = plasma_get_numthreads();


    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zheev_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zheev_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zheev_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zheev_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(&descT) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zheev_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /*
    if (jobz == PlasmaVec){
       if (plasma_desc_check(Q) != PLASMA_SUCCESS) {
           plasma_error("PLASMA_zheev_Tile_Async", "invalid descriptor");
           return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
       }
    }*/
    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_zheev_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_zheev_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        plasma_error("PLASMA_zheev_Tile_Async", "matrix need to be square");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zheev_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /*
    if (jobz == PlasmaVec) {
        plasma_error("PLASMA_zheev_Tile_Async", "computing the eigenvectors is not supported in this version");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if (jobz == PlasmaVec){
       if (Q->nb != Q->mb) {
           plasma_error("PLASMA_zheev_Tile_Async", "only square tiles supported");
           return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
       }
    }
   */
    AB = (PLASMA_Complex64_t *)plasma_shared_alloc(plasma, LDAB*descA.ln, PlasmaComplexDouble);
    if (AB == NULL) {
        plasma_error("PLASMA_zheev", "plasma_shared_alloc() failed");
        plasma_shared_free(plasma, AB);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    E = (double *)plasma_shared_alloc(plasma, descA.n-1, PlasmaRealDouble);
    if (E == NULL) {
        plasma_error("PLASMA_zheev", "plasma_shared_alloc() failed");
        plasma_shared_free(plasma, E);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }


    /* Currently NOT equivalent to LAPACK's
     */

    /* Reduction to tridiagonal form
     * with a two-stage approach.
     */

    double timelpk=0.0,tconvert=0.0,timeaplQ2=0.0,timeaplQ=0.0, timeeigen=0.0, timeT=0.0;
    double timeB=0.0,timeblg=0.0,timeaplQ1=0.0,timeconv1=0.0,timeconv2=0.0,timeall=0.0;
    timeall = PLASMA_Wtime();
    timeB   = PLASMA_Wtime();

    /*=======================================
     *  calling Reduction to BAND
     *  then convert matrix to band form
     *=======================================*/
    plasma_dynamic_call_5(plasma_pzhetrd_he2hb,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_6( plasma_pzhbcpy_t2bl,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_Complex64_t, AB,
        int, LDAB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_sync();

    stat = sequence->status;
    if (stat != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zheev","pzhetrd_he2hb+pzcopy");
        return stat;
    }

    timeB   = PLASMA_Wtime()-timeB;
    printf("\n  Finish reduction to band time= %lf \n",timeB);
    /*=======================================
     *  END of calling Reduction to BAND
     *=======================================*/


    /*=======================================
     *  calling bulge chasing
     *=======================================*/
    PLASMA_Complex64_t *TAU2;
    PLASMA_Complex64_t *V2;
    PLASMA_Complex64_t *T2;
    int Vblksiz,blkcnt,LDT,LDV,NBTILES, mklth;
    int WANTZ   = 0;
    int blguplo = PlasmaLower;
    int NE      = N;
    if(jobz==PlasmaNoVec)
            WANTZ=0;
    else
            WANTZ=2;

    Vblksiz  = min(NB,40);
    LDT      = Vblksiz;
    NBTILES  = descA.lnt;
    NBTILES  = plasma_ceildiv(N,NB);

    if(jobz == PlasmaVec){
        findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
        TAU2   = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, blkcnt*Vblksiz,     PlasmaComplexDouble);
        V2     = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, blkcnt*LDV*Vblksiz, PlasmaComplexDouble);
        T2     = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, blkcnt*LDT*Vblksiz, PlasmaComplexDouble);
        if ( (TAU2 == NULL) || (V2 == NULL) || (T2 == NULL) ) {
            plasma_error("PLASMA_zheev", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            plasma_shared_free(plasma, T2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0, blkcnt*Vblksiz*sizeof(PLASMA_Complex64_t));
        memset(V2,   0, blkcnt*LDV*Vblksiz*sizeof(PLASMA_Complex64_t));
        memset(T2,   0, blkcnt*LDT*Vblksiz*sizeof(PLASMA_Complex64_t));
    }else{
        TAU2   = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*N, PlasmaComplexDouble);
        V2     = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*N, PlasmaComplexDouble);
        if ( (TAU2 == NULL) || (V2 == NULL) ) {
            plasma_error("PLASMA_zheev", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0, 2*N*sizeof(PLASMA_Complex64_t));
        memset(V2,   0, 2*N*sizeof(PLASMA_Complex64_t));
    }


    timeblg  = PLASMA_Wtime();
    plasma_static_call_13(plasma_pzhetrd_hb2st_v1,
        PLASMA_enum,         blguplo,
        int,                 N,
        int,                 NB,
        int,                 Vblksiz,
        PLASMA_Complex64_t*, AB,
        int,                 LDAB,
        PLASMA_Complex64_t*, V2,
        PLASMA_Complex64_t*, TAU2,
        double*,             W,
        double*,             E,
        int,                 WANTZ,
        PLASMA_sequence*,    sequence,
        PLASMA_request*,     request);

    plasma_dynamic_sync();
    timeblg  = PLASMA_Wtime()-timeblg;
    printf("  Finish Bulge       timing= %lf \n", timeblg);
    /*=======================================
     *  END of calling bulge chasing
     *=======================================*/

    /*=======================================
     *  calling eigensolver
     *=======================================*/
    plasma_unsetaffinity();
    if(jobz == PlasmaVec){
       memset(Q,0,N*LDQ*sizeof(PLASMA_Complex64_t));
    }
    printf("  start eigensolver\n");
    timeeigen = PLASMA_Wtime();
    timelpk = PLASMA_Wtime();
#if defined(USEMKL)
    mklth=THREADS; // min(12,THREADS)
    mkl_set_num_threads(mklth);
#endif
    // call eigensolver using lapack routine for our resulting tridiag [D E]
    if(jobz == PlasmaNoVec){
        core_zstedc(PlasmaNoVec, N, W, E, Q, LDQ);
    }else{
        core_zstedc(PlasmaIvec, N, W, E, Q, LDQ);
    }
#if defined(USEMKL)
    mkl_set_num_threads( 1 );
#endif
    timelpk = PLASMA_Wtime()-timelpk;
    printf("  Finish Eigensolver timing= %lf  WANTZ %d threads %d\n", timelpk, WANTZ, mklth);
    plasma_setaffinity(0);
    /*=======================================
     *  END of calling eigensolver
     *=======================================*/









/*
    FILE *trace_file;
    trace_file = fopen("AJETE/Q2", "w");
    for (j = 0; j < N ; j++)
          for (i = 0; i < N ; i++)
                         fprintf(trace_file,"%10d%10d%40.30e\n",i+1,j+1,Q[j*N+i]);
    fclose(trace_file);
*/



     if (jobz == PlasmaVec){
        /*=======================================
         *  apply Q2 from the bulge
         *=======================================*/
        // compute in parallel T2 from bulge
        timeT  = PLASMA_Wtime();
        plasma_static_call_8(plasma_pzlarft_blgtrd,
            int,                 N,
            int,                 NB,
            int,                 Vblksiz,
            PLASMA_Complex64_t*, V2,
            PLASMA_Complex64_t*, T2,
            PLASMA_Complex64_t*, TAU2,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        plasma_dynamic_sync();
        timeT  = PLASMA_Wtime()-timeT;
        printf("  Finish compute T2  timing= %lf \n", timeT);


        timeaplQ2  = PLASMA_Wtime();
        plasma_static_call_14(plasma_pzunmqr_blgtrd,
            PLASMA_enum,         PlasmaLeft,
            PLASMA_enum,         PlasmaNoTrans,
            int,                 N,
            int,                 NB,
            int,                 N,
            int,                 Vblksiz,
            int,                 WANTZ,
            PLASMA_Complex64_t*, V2,
            PLASMA_Complex64_t*, T2,
            PLASMA_Complex64_t*, TAU2,
            PLASMA_Complex64_t*, Q,
            int,                 LDQ,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        plasma_dynamic_sync();
        timeaplQ2  = PLASMA_Wtime()-timeaplQ2;
        printf("  Finish compute Q2  timing= %lf \n", timeaplQ2);





        /*=======================================
         *  apply Q1 from the reduction to band
         *=======================================*/
         /*CASE NB>N no Q1 to be applied only bulge chasing has been done*/
         if(NB<N){
         plasma_dynamic_sync();
         timeconv1   = PLASMA_Wtime();
         plasma_zooplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, N, N, plasma_desc_mat_free(&(descQ)) );
        // plasma_ziplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, N, N );
         plasma_dynamic_sync();
         timeconv1    = PLASMA_Wtime()-timeconv1;
         timeaplQ1   = PLASMA_Wtime();
         // Accumulate the transformations from the first stage
         //plasma_parallel_call_7(plasma_pzunmqr, // static code is not available for Notrans
         if(uplo==PlasmaLower){
             plasma_dynamic_call_7(plasma_pzunmqr,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaNoTrans,
                PLASMA_desc, plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb),
                PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n),
                PLASMA_desc, plasma_desc_submatrix(descT, descT.mb, 0, descT.m-descT.mb, descT.n-descT.nb),
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);

         }else if(uplo==PlasmaUpper){
             plasma_dynamic_call_7(plasma_pzunmlq,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaConjTrans,
                PLASMA_desc, plasma_desc_submatrix(descA, 0, descA.nb, descA.n-descA.nb, descA.m-descA.mb),
                PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n         ),
                PLASMA_desc, plasma_desc_submatrix(descT, 0, descT.nb, descT.n-descT.nb, descT.m-descT.mb),
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
         }
         plasma_dynamic_sync();
         timeaplQ1   = PLASMA_Wtime()-timeaplQ1;
         printf("  Finish compute Q1  timing= %lf \n", timeaplQ1);
         timeconv2    = PLASMA_Wtime();
         plasma_zooptile2lap( descQ, Q, NB, NB, LDQ, N );
         //plasma_ziptile2lap( descQ, Q, NB, NB, LDQ, N );
         plasma_dynamic_sync();
         timeconv2    = PLASMA_Wtime()-timeconv2;
         printf("  Finish convert     timing= %lf \n", timeconv1+timeconv2);

         /*
         plasma_dynamic_sync();
         descQ = plasma_desc_init( PlasmaComplexDouble, NB, NB, NB*NB, LDQ, N, 0, 0, N, N);
         descQ = plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n);
         descQ.mat = Q+NB;
         timeaplQ1   = PLASMA_Wtime();
         // Accumulate the transformations from the first stage
         if(uplo==PlasmaLower){
             plasma_dynamic_call_7(plasma_pzunmqr_tlpk,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaNoTrans,
                PLASMA_desc, plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb),
                PLASMA_desc, descQ,
                PLASMA_desc, plasma_desc_submatrix(descT, descT.mb, 0, descT.m-descT.mb, descT.n-descT.nb),
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);

         }else if(uplo==PlasmaUpper){
             plasma_static_call_7(plasma_pzunmlq_tlpk,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaConjTrans,
                PLASMA_desc, plasma_desc_submatrix(descA, 0, descA.nb, descA.n-descA.nb, descA.m-descA.mb),
                PLASMA_desc, descQ,
                PLASMA_desc, plasma_desc_submatrix(descT, 0, descT.nb, descT.n-descT.nb, descT.m-descT.mb),
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
         }
         plasma_dynamic_sync();
         timeaplQ1   = PLASMA_Wtime()-timeaplQ1;
         printf("  Finish compute Q1  timing= %lf \n", timeaplQ1);
         printf("  Finish convert     timing= %lf \n", timeconv1+timeconv2);

         */
         //======================
         }// END of if(NB<N)
         //======================

        /* Set the V's to zero before the 2nd stage (bulge chasing) */
     /*
    int NN=descA.lm;
    PLASMA_Complex64_t *QQ=descQ.mat;
    FILE *trace_file;
    printf("voici NN %d\n",NN);
    trace_file = fopen("AJETE/Q", "w");
    for (j = 0; j < NN ; j++)
          for (i = 0; i < NN ; i++)
                         fprintf(trace_file,"%10d%10d%40.30e\n",i+1,j+1,Q[j*NN+i]);
    fclose(trace_file);

    trace_file = fopen("AJETE/D", "w");
    for (j = 0; j < NN ; j++)
                         fprintf(trace_file,"%10d%10d%40.30e\n",j+1,1,W[j]);
    fclose(trace_file);
*/

     }


    timeall = PLASMA_Wtime()-timeall;
    printf("  Finish full eigenproblem threads %d  N %d  timeall= %lf \n",THREADS, N, timeall);


    plasma_shared_free(plasma, E);

    return PLASMA_SUCCESS;
}
