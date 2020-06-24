/**
 * @file zhetrd.c
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
#include "common.h"
/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zhetrd - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary
 *  Q2).  Let Q = Q1 * Q2 be the global unitary transformation; Q**H *
 *  A * Q = S.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes tridiagonal only;
 *          = PlasmaVec: computes tridiagonal and generate the orthogonal matrix Q.
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
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix:
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix:
 *          E(i) = A(i,i+1) if uplo = PlasmaUpper, E(i) = A(i+1,i) if uplo = PlasmaLower.
 *
 * @param[out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_zhetrd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec, then if return value = 0, Q
 *          contains the N-by-N unitary matrix Q.
 *          If jobz = PlasmaNoVec, then it is not referenced.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= N.
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
 * @sa PLASMA_zhetrd_Tile
 * @sa PLASMA_zhetrd_Tile_Async
 * @sa PLASMA_chetrd
 * @sa PLASMA_dsytrd
 * @sa PLASMA_ssytrd
 *
 ******************************************************************************/
int PLASMA_zhetrd(PLASMA_enum jobz, PLASMA_enum uplo, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 double *D,
                 double *E,
                 PLASMA_desc *descT,
                 PLASMA_Complex64_t *Q, int LDQ)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_zhetrd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_zhetrd", "illegal value of jobz");
        return -1;
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_zhetrd", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zhetrd", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zhetrd", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZHETRD, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhetrd", "plasma_tune() failed");
        return status;
    }
    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zhetrd_Tile_Async(jobz, uplo, &descA, D, E, descT, Q, LDQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zhetrd_Tile - reduces a complex Hermitian matrix A to real symmetric
 *  tridiagonal form S using a two-stage approach
 *  First stage: reduction to band tridiagonal form (unitary Q1);
 *  Second stage: reduction from band to tridiagonal form (unitary Q2).
 *  Let Q = Q1 * Q2 be the global unitary transformation;
 *  Q**H * A * Q = S.
 *  Tile equivalent of PLASMA_zhetrd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes tridiagonal only;
 *          = PlasmaVec: computes tridiagonal and generate the orthogonal matrix Q.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.  If uplo
 *          = PlasmaUpper, the leading N-by-N upper triangular part of
 *          A contains the upper triangular part of the matrix A, and
 *          the strictly lower triangular part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] D
 *          On exit, the diagonal elements of the tridiagonal matrix:
 *          D(i) = A(i,i).
 *
 * @param[out] E
 *          On exit, he off-diagonal elements of the tridiagonal matrix:
 *          E(i) = A(i,i+1) if uplo = PlasmaUpper,
 *          E(i) = A(i+1,i) if uplo = PlasmaLower.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec, then if return value = 0, Q
 *          contains the N-by-N unitary matrix Q.
 *          If jobz = PlasmaNoVec, then it is not referenced.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= N.
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
 * @sa PLASMA_zhetrd
 * @sa PLASMA_zhetrd_Tile_Async
 * @sa PLASMA_chetrd_Tile
 * @sa PLASMA_dsytrd_Tile
 * @sa PLASMA_ssytrd_Tile
 * @sa PLASMA_zhetrd_Tile
 *
 ******************************************************************************/
int PLASMA_zhetrd_Tile(PLASMA_enum jobz, PLASMA_enum uplo,
                       PLASMA_desc *A, double *D, double *E,
                       PLASMA_desc *T, PLASMA_Complex64_t *Q, int LDQ)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zhetrd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zhetrd_Tile_Async(jobz, uplo, A, D, E, T, Q, LDQ, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zhetrd_Tile_Async - Computes all eigenvalues and,
 *  optionally, eigenvectors of a complex Hermitian matrix A using a
 *  two-stage approach:
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
 * @sa PLASMA_zhetrd
 * @sa PLASMA_zhetrd_Tile
 * @sa PLASMA_chetrd_Tile_Async
 * @sa PLASMA_dsytrd_Tile_Async
 * @sa PLASMA_ssytrd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zhetrd_Tile_Async(PLASMA_enum jobz,
                            PLASMA_enum uplo,
                            PLASMA_desc *A,
                            double *W,
                            double *E,
                            PLASMA_desc *T,
                            PLASMA_Complex64_t *Q, int LDQ,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descQ;
    PLASMA_Complex64_t *AB;
    int N, NB, LDAB;
    int i;
    int status;


    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zhetrd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zhetrd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zhetrd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhetrd_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhetrd_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_zhetrd_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_zhetrd_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        plasma_error("PLASMA_zhetrd_Tile_Async", "matrix need to be square");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zhetrd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    N    = descA.m;
    NB   = min(descA.mb,descA.m);
    LDAB = 2*NB+1;

    /* Allocate workspace for band storage of the band matrix A and for the off diagonal after tridiagonalisation */
    AB = (PLASMA_Complex64_t *)plasma_shared_alloc(plasma, LDAB*N, PlasmaComplexDouble);
    memset( AB, 0, LDAB * N * sizeof(PLASMA_Complex64_t) );
    if (AB == NULL) {
        plasma_error("PLASMA_zhetrd", "plasma_shared_alloc() failed");
        plasma_shared_free(plasma, AB);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    /* Reduction to tridiagonal form
     * with a two-stage approach.
     */
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
        PLASMA_Complex64_t*, AB,
        int, LDAB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_sync();

    status = sequence->status;
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhetrd","pzhetrd_he2hb+pzcopy");
        return status;
    }
    /*=======================================
     *  END of calling Reduction to BAND
     *=======================================*/
    /*=======================================
     *  calling bulge chasing
     *=======================================*/
    PLASMA_Complex64_t *TAU2 = NULL;
    PLASMA_Complex64_t *V2 = NULL;
    PLASMA_Complex64_t *T2 = NULL;
    int Vblksiz, blkcnt, LDT, LDV;
    int WANTZ   = 0;
    int blguplo = PlasmaLower;
    /* int NE      = N; // for later use when a portion of the eigenvectors are requested*/
    if( jobz == PlasmaNoVec )
        WANTZ=0;
    else
        WANTZ=1;

    /* Vblksiz correspond to the blocking used when applying V2 to the matrix Q
     * it is similar to IB in LAPACK ZUNMQR.
     * blkcnt is the number of losange or tile of Vs */
    /* Note that in case PlamaVec requested, the V2 and T2 are stored by the
     * bulgechasing function in a special format:
     * for V2s: it store the V2(LDV,Vblksiz) of each losange in a tile storage meaning
     * that V2_1 is stored then V2_2,..., V2_blkcnt.
     * blkcnt is the number of losange.
     * */
    Vblksiz = min(NB,48);
    LDT     = Vblksiz;
    if( jobz == PlasmaVec ) {
        findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
        TAU2   = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma,     blkcnt*Vblksiz, PlasmaComplexDouble);
        V2     = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, LDV*blkcnt*Vblksiz, PlasmaComplexDouble);
        T2     = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, LDT*blkcnt*Vblksiz, PlasmaComplexDouble);
        if ( (TAU2 == NULL) || (V2 == NULL) || (T2 == NULL) ) {
            plasma_error("PLASMA_zhetrd", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            plasma_shared_free(plasma, T2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0,     blkcnt*Vblksiz*sizeof(PLASMA_Complex64_t));
        memset(V2,   0, LDV*blkcnt*Vblksiz*sizeof(PLASMA_Complex64_t));
        memset(T2,   0, LDT*blkcnt*Vblksiz*sizeof(PLASMA_Complex64_t));
    }
    else {
        TAU2   = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*N, PlasmaComplexDouble);
        V2     = (PLASMA_Complex64_t *) plasma_shared_alloc(plasma, 2*N, PlasmaComplexDouble);
        if ( (TAU2 == NULL) || (V2 == NULL) ) {
            plasma_error("PLASMA_zhetrd", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0, 2*N*sizeof(PLASMA_Complex64_t));
        memset(V2,   0, 2*N*sizeof(PLASMA_Complex64_t));
    }

    plasma_parallel_call_13(plasma_pzhetrd_hb2st_v1,
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
    /* WARNING: If plasma_pzhetrd_hb2st is implemented through a dynamic call, don't
     * forget to synchronize */
    plasma_dynamic_sync();
    /*=======================================
     *  END of calling bulge chasing
     *=======================================*/
    if (jobz == PlasmaVec){
        /*=======================================
         *  generate Q2 from the bulge
         *=======================================*/
        /* Initialize Q to Identity */
        memset(Q,   0, N*LDQ*sizeof(PLASMA_Complex64_t));
        for(i=0; i<N; i++){
            Q[i+i*LDQ] = 1.0;
        }
        /* compute T2 */
        plasma_static_call_8(plasma_pzlarft_blgtrd,
            int,                 N,
            int,                 NB,
            int,                 Vblksiz,
            PLASMA_Complex64_t*, V2,
            PLASMA_Complex64_t*, T2,
            PLASMA_Complex64_t*, TAU2,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);
        /* apply Q2 from Left */
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
        /*=======================================
         *  apply Q1 from the reduction to band
         *=======================================*/
        /* CASE NB>N, Q1 doesn't need to be applied, only bulge chasing has been done */
        if( NB < N ){
            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_zooplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, N, N, sequence, request, plasma_desc_mat_free(&(descQ)) );
            } else {
                plasma_ziplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, N, N, sequence, request );
            }

            /* Accumulate the transformations from the first stage */
            if(uplo==PlasmaLower){
                plasma_dynamic_call_7(plasma_pzunmqr,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n),
                    PLASMA_desc, plasma_desc_submatrix(descT, descT.mb, 0, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);

            }
            else {
                plasma_dynamic_call_7(plasma_pzunmlq,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaConjTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, 0, descA.nb, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n         ),
                    PLASMA_desc, plasma_desc_submatrix(descT, 0, descT.nb, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_zooptile2lap( descQ, Q, NB, NB, LDQ, N, sequence, request );
                plasma_dynamic_sync();
                plasma_desc_mat_free(&descQ);
            } else {
                plasma_ziptile2lap( descQ, Q, NB, NB, LDQ, N, sequence, request );
                plasma_dynamic_sync();
            }
        } /* END of ( NB < N ) */
    }
    /*=======================================
     *  END of calling computing Q
     *=======================================*/
    if( jobz == PlasmaVec )
        plasma_shared_free(plasma, T2);
    plasma_shared_free(plasma, V2);
    plasma_shared_free(plasma, TAU2);
    plasma_shared_free(plasma, AB);
    return PLASMA_SUCCESS;
}
