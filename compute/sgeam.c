/**
 *
 * @file sgeam.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_sgeam - Performs one of the matrix-matrix operations
 *
 *    \f[ B = \alpha [op( A ) + \beta opt( B ) \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = g( X' )
 *
 *  alpha and beta are scalars, and A, B  are matrices, with op( A )
 *  an m by n matrix, and opt(B) an m by n matrix.
 *  B = A + B; B = conj(A') + B; conj(B') = A + conj(B'); conj(B') = conj(A') + conj (B');
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( A ) and of the matrix C. N >= 0.
 *
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgeam_Tile
 * @sa PLASMA_cgeam
 * @sa PLASMA_sgeam
 * @sa PLASMA_sgeam
 *
 ******************************************************************************/
int PLASMA_sgeam(PLASMA_enum transA, PLASMA_enum transB, int M, int N,
                 float alpha, float *A, int LDA, float beta,  float *B, int LDB)
{
    int NB;
    int Am, An, Bm, Bn;
    int status;
    PLASMA_desc descA, descB;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgeam", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((transA != PlasmaNoTrans) && (transA != PlasmaTrans) && (transA != PlasmaConjTrans)) {
        plasma_error("PLASMA_sgeam", "illegal value of transA");
        return -1;
    }
    if ((transB != PlasmaNoTrans) && (transB != PlasmaTrans) && (transB != PlasmaConjTrans)) {
        plasma_error("PLASMA_sgeam", "illegal value of transB");
        return -2;
    }
    if ( transA == PlasmaNoTrans ) { 
        Am = M; An = N;
    } else {
        Am = N; An = M;
    }
    if ( transB == PlasmaNoTrans ) { 
        Bm = M; Bn = N;
    } else {
        Bm = N; Bn = M;
    }
    if (M < 0) {
        plasma_error("PLASMA_sgeam", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_sgeam", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_sgeam", "illegal value of LDA");
        return -8;
    }
    if (LDB < max(1, Bm)) {
        plasma_error("PLASMA_sgeam", "illegal value of LDB");
        return -10;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (float)0.0) && beta == (float)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGEAM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeam", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_sooplap2tile( descB, B, NB, NB, LDB, Bn, 0, 0, Bm, Bn, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, 
                            sequence, &request);
        plasma_siplap2tile( descB, B, NB, NB, LDB, Bn, 0, 0, Bm, Bn, 
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_sgeam_Tile_Async(
        transA, transB, alpha, &descA, beta, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, An,  sequence, &request);
        plasma_siptile2lap( descB, B, NB, NB, LDB, Bn,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_sgeam_Tile - Performs matrix summation.
 *  Tile equivalent of PLASMA_sgeam().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-na matrix, where na is N when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a LDB-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A ) + beta*opt( B ) )
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgeam
 * @sa PLASMA_sgeam_Tile_Async
 * @sa PLASMA_cgeam_Tile
 * @sa PLASMA_sgeam_Tile
 * @sa PLASMA_sgeam_Tile
 *
 ******************************************************************************/
int PLASMA_sgeam_Tile(PLASMA_enum transA, PLASMA_enum transB,
                      float alpha, PLASMA_desc *A, float beta, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;
    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgemm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sgeam_Tile_Async(transA, transB, alpha, A, beta, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sgeam_Tile_Async - Performs matrix summation.
 *  Non-blocking equivalent of PLASMA_sgeam_Tile().
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
 * @sa PLASMA_sgeam
 * @sa PLASMA_sgeam_Tile
 * @sa PLASMA_cgeam_Tile_Async
 * @sa PLASMA_sgeam_Tile_Async
 * @sa PLASMA_sgeam_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sgeam_Tile_Async(PLASMA_enum transA, PLASMA_enum transB,
                            float alpha, PLASMA_desc *A, float beta, PLASMA_desc *B,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descB;
 
    int M, N;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;
    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgeam_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sgeam_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sgeam_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeam_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgeam_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if ((transA != PlasmaNoTrans) && (transA != PlasmaTrans) && (transA != PlasmaConjTrans)) {
        plasma_error("PLASMA_sgeam_Tile_Async", "illegal value of transA");
        return plasma_request_fail(sequence, request, -1);
    }
    if ((transB != PlasmaNoTrans) && (transB != PlasmaTrans) && (transB != PlasmaConjTrans)) {
        plasma_error("PLASMA_sgeam_Tile_Async", "illegal value of transB");
        return plasma_request_fail(sequence, request, -2);
    }

    if ( transA == PlasmaNoTrans ) {
        Am  = descA.m;
        An  = descA.n;
        Amb = descA.mb;
        Anb = descA.nb;
        Ai  = descA.i;
        Aj  = descA.j;
    } else {
        Am  = descA.n;
        An  = descA.m;
        Amb = descA.nb;
        Anb = descA.mb;
        Ai  = descA.j;
        Aj  = descA.i;
    }

    if ( transB == PlasmaNoTrans ) {
        Bm  = descB.m;
        Bn  = descB.n;
        Bmb = descB.mb;
        Bnb = descB.nb;
        Bi  = descB.i;
        Bj  = descB.j;
    } else {
        Bm  = descB.n;
        Bn  = descB.m;
        Bmb = descB.nb;
        Bnb = descB.mb;
        Bi  = descB.j;
        Bj  = descB.i;
    }

    if ( (Amb != descB.mb) || (Anb != descB.nb) ) {
        plasma_error("PLASMA_sgeam_Tile_Async", "tile sizes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai != descB.i) || (Aj != Bi) || (Bj != descA.j) ) {
        plasma_error("PLASMA_sgeam_Tile_Async", "start indexes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    M = descA.m;
    N = descA.n;
    

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (float)0.0 ) && beta == (float)1.0))
        return PLASMA_SUCCESS;
    plasma_parallel_call_9(plasma_psgeam,
        PLASMA_enum, transA,
        PLASMA_enum, transB,
        float, alpha,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        float, beta,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
