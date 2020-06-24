/**
 *
 * @file coreblas_z.c
 *
 *  PLASMA core_blas tracing kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This file provides the wrapper for each function of the
 *  core_blas library which will generate an event before and
 *  after the execution of the kernel.
 *  This file is automatically generated with convert2eztrace.pl
 *  script. DO NOT MANUALLY EDIT THIS FILE.
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <eztrace.h>
#include <ev_codes.h>
#include "common.h"
#include "coreblas_ev_codes.h"
#include "coreblas_macros.h"
#undef REAL
#define COMPLEX
FUNCTION_VOID( CORE_dzasum, ASUM, void ,
          (PLASMA_enum storev, PLASMA_enum uplo, int M, int N, const PLASMA_Complex64_t *A, int lda, double *work),
          (storev, uplo, M, N, A, lda, work) )
FUNCTION_VOID( CORE_zbrdalg1, BRDALG, void ,
          ( PLASMA_enum uplo, int n, int nb, PLASMA_Complex64_t *A, int lda, PLASMA_Complex64_t *VQ, PLASMA_Complex64_t *TAUQ, PLASMA_Complex64_t *VP, PLASMA_Complex64_t *TAUP, int Vblksiz, int wantz, int i, int sweepid, int m, int grsiz, PLASMA_Complex64_t *work),
          (uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP, Vblksiz, wantz, i, sweepid, m, grsiz, work) )
FUNCTION_TYPE( CORE_zgeadd, GEADD, int ,
          (int M, int N, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB),
          (M, N, alpha, A, LDA, B, LDB) )
FUNCTION_TYPE( CORE_zgelqt, GELQT, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK),
          (M, N, IB, A, LDA, T, LDT, TAU, WORK) )
FUNCTION_VOID( CORE_zgemm, GEMM, void ,
          (PLASMA_enum transA, int transB, int M, int N, int K, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC),
          (transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) )
FUNCTION_VOID( CORE_zgemv, GEMV, void ,
          (PLASMA_enum trans, int m, int n, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda, const PLASMA_Complex64_t *x, int incx, PLASMA_Complex64_t beta, PLASMA_Complex64_t *y, int incy),
          (trans, m, n, alpha, A, lda, x, incx, beta, y, incy) )
FUNCTION_VOID( CORE_zgeqp3_init, GEQP3_INIT, void ,
          ( int n, int *jpvt ),
          (n, jpvt) )
FUNCTION_VOID( CORE_zgeqp3_larfg, LARFG, void ,
          ( PLASMA_desc A, int ii, int jj, int i, int j, PLASMA_Complex64_t *tau, PLASMA_Complex64_t *beta ),
          (A, ii, jj, i, j, tau, beta) )
FUNCTION_VOID( CORE_zgeqp3_norms, GEQP3_NORMS, void ,
          ( PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 ),
          (A, ioff, joff, norms1, norms2) )
FUNCTION_VOID( CORE_zgeqp3_pivot, GEQP3_PIVOT, void ,
          ( PLASMA_desc A, PLASMA_Complex64_t *F, int ldf, int jj, int k, int *jpvt, double *norms1, double *norms2, int *info ),
          (A, F, ldf, jj, k, jpvt, norms1, norms2, info) )
FUNCTION_TYPE( CORE_zgeqp3_tntpiv, GEQRT, int ,
          (int m, int n, PLASMA_Complex64_t *A, int lda, int *IPIV, PLASMA_Complex64_t *tau, int *iwork),
          (m, n, A, lda, IPIV, tau, iwork) )
FUNCTION_VOID( CORE_zgeqp3_update, GEQP3_UPDATE, void ,
          ( const PLASMA_Complex64_t *Ajj, int lda1, PLASMA_Complex64_t *Ajk, int lda2, const PLASMA_Complex64_t *Fk, int ldf, int joff, int k, int koff, int nb, double *norms1, double *norms2, int *info ),
          (Ajj, lda1, Ajk, lda2, Fk, ldf, joff, k, koff, nb, norms1, norms2, info) )
FUNCTION_TYPE( CORE_zgeqrt, GEQRT, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK),
          (M, N, IB, A, LDA, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_zgessm, GESSM, int ,
          (int M, int N, int K, int IB, const int *IPIV, const PLASMA_Complex64_t *L, int LDL, PLASMA_Complex64_t *A, int LDA),
          (M, N, K, IB, IPIV, L, LDL, A, LDA) )
FUNCTION_TYPE( CORE_zgessq, LASSQ, int ,
          (int M, int N, const PLASMA_Complex64_t *A, int LDA, double *scale, double *sumsq),
          (M, N, A, LDA, scale, sumsq) )
FUNCTION_TYPE( CORE_zgetrf, GETRF, int ,
          (int m, int n, PLASMA_Complex64_t *A, int lda, int *IPIV, int *info),
          (m, n, A, lda, IPIV, info) )
FUNCTION_TYPE( CORE_zgetrf_incpiv, GETRF, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A, int LDA, int *IPIV, int *INFO),
          (M, N, IB, A, LDA, IPIV, INFO) )
FUNCTION_TYPE( CORE_zgetrf_nopiv, GETRF, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A, int LDA),
          (M, N, IB, A, LDA) )
FUNCTION_TYPE( CORE_zgetrf_reclap, GETRF, int ,
          (int M, int N, PLASMA_Complex64_t *A, int LDA, int *IPIV, int *info),
          (M, N, A, LDA, IPIV, info) )
FUNCTION_TYPE( CORE_zgetrf_rectil, GETRF, int ,
          (const PLASMA_desc A, int *IPIV, int *info),
          (A, IPIV, info) )
FUNCTION_VOID( CORE_zgetrip, GETRIP, void ,
          (int m, int n, PLASMA_Complex64_t *A, PLASMA_Complex64_t *W) ,
          (m, n, A, W)  )
FUNCTION_VOID( CORE_zhegst, HEGST, void ,
          (int itype, PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB, int *INFO),
          (itype, uplo, N, A, LDA, B, LDB, INFO) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_zhemm, HEMM, void ,
          (PLASMA_enum side, PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC),
          (side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC) )
#endif
#ifdef COMPLEX
FUNCTION_VOID( CORE_zher2k, HER2K, void ,
          (PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *B, int LDB, double beta, PLASMA_Complex64_t *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) )
#endif
FUNCTION_TYPE( CORE_zherfb, HERFB, int ,
          ( PLASMA_enum uplo, int n, int k, int ib, int nb, const PLASMA_Complex64_t *A, int lda, const PLASMA_Complex64_t *T, int ldt, PLASMA_Complex64_t *C, int ldc, PLASMA_Complex64_t *WORK, int ldwork ),
          (uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_zherk, HERK, void ,
          (PLASMA_enum uplo, PLASMA_enum trans, int N, int K, double alpha, const PLASMA_Complex64_t *A, int LDA, double beta, PLASMA_Complex64_t *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, beta, C, LDC) )
#endif
#ifdef COMPLEX
FUNCTION_TYPE( CORE_zhessq, LASSQ, int ,
          (PLASMA_enum uplo, int N, const PLASMA_Complex64_t *A, int LDA, double *scale, double *sumsq),
          (uplo, N, A, LDA, scale, sumsq) )
#endif
FUNCTION_VOID( CORE_zlacpy, LACPY, void ,
          (PLASMA_enum uplo, int M, int N, const PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB),
          (uplo, M, N, A, LDA, B, LDB) )
FUNCTION_TYPE( CORE_zlacpy_pivot, LACPY, int ,
          ( const PLASMA_desc descA, PLASMA_enum direct, int k1, int k2, const int *ipiv, int *rankin, int *rankout, PLASMA_Complex64_t *A, int lda, int init),
          (descA, direct, k1, k2, ipiv, rankin, rankout, A, lda, init) )
FUNCTION_VOID( CORE_zlange, LANGE, void ,
          (int norm, int M, int N, const PLASMA_Complex64_t *A, int LDA, double *work, double *normA),
          (norm, M, N, A, LDA, work, normA) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_zlanhe, LANHE, void ,
          (int norm, PLASMA_enum uplo, int N, const PLASMA_Complex64_t *A, int LDA, double *work, double *normA),
          (norm, uplo, N, A, LDA, work, normA) )
#endif
FUNCTION_VOID( CORE_zlansy, LANSY, void ,
          (int norm, PLASMA_enum uplo, int N, const PLASMA_Complex64_t *A, int LDA, double *work, double *normA),
          (norm, uplo, N, A, LDA, work, normA) )
FUNCTION_VOID( CORE_zlantr, LANGE, void ,
          (PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N, const PLASMA_Complex64_t *A, int LDA, double *work, double *normA),
          (norm, uplo, diag, M, N, A, LDA, work, normA) )
FUNCTION_TYPE( CORE_zlarfb_gemm, LARFB, int ,
          (PLASMA_enum side, PLASMA_enum trans, int direct, int storev, int M, int N, int K, const PLASMA_Complex64_t *V, int LDV, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *C, int LDC, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, direct, storev, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK) )
FUNCTION_VOID( CORE_zlaset2, LASET, void ,
          (PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA),
          (uplo, M, N, alpha, A, LDA) )
FUNCTION_VOID( CORE_zlaset, LASET, void ,
          (PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta, PLASMA_Complex64_t *A, int LDA),
          (uplo, M, N, alpha, beta, A, LDA) )
FUNCTION_VOID( CORE_zlaswp, LASWP, void ,
          (int N, PLASMA_Complex64_t *A, int LDA, int I1, int I2, const int *IPIV, int INC),
          (N, A, LDA, I1, I2, IPIV, INC) )
FUNCTION_TYPE( CORE_zlaswp_ontile, LASWP, int ,
          (PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc),
          (descA, i1, i2, ipiv, inc) )
FUNCTION_TYPE( CORE_zswptr_ontile, TRSM, int ,
          (PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc, const PLASMA_Complex64_t *Akk, int ldak),
          (descA, i1, i2, ipiv, inc, Akk, ldak) )
FUNCTION_TYPE( CORE_zlaswpc_ontile, LASWP, int ,
          (PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc),
          (descA, i1, i2, ipiv, inc) )
FUNCTION_TYPE( CORE_zlatro, LATRO, int ,
          (PLASMA_enum uplo, PLASMA_enum trans, int M, int N, const PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB),
          (uplo, trans, M, N, A, LDA, B, LDB) )
FUNCTION_VOID( CORE_zlauum, LAUUM, void ,
          (PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA),
          (uplo, N, A, LDA) )
FUNCTION_TYPE( CORE_zpemv, PEMV, int ,
          (PLASMA_enum trans, int storev, int M, int N, int L, PLASMA_Complex64_t ALPHA, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *X, int INCX, PLASMA_Complex64_t BETA, PLASMA_Complex64_t *Y, int INCY, PLASMA_Complex64_t *WORK),
          (trans, storev, M, N, L, ALPHA, A, LDA, X, INCX, BETA, Y, INCY, WORK) )
#ifdef COMPLEX
FUNCTION_VOID( CORE_zplghe, PLGHE, void ,
          ( double bump, int m, int n, PLASMA_Complex64_t *A, int lda, int bigM, int m0, int n0, unsigned long long int seed ),
          (bump, m, n, A, lda, bigM, m0, n0, seed) )
#endif
FUNCTION_VOID( CORE_zplgsy, PLGSY, void ,
          ( PLASMA_Complex64_t bump, int m, int n, PLASMA_Complex64_t *A, int lda, int bigM, int m0, int n0, unsigned long long int seed ),
          (bump, m, n, A, lda, bigM, m0, n0, seed) )
FUNCTION_VOID( CORE_zplrnt, PLRNT, void ,
          ( int m, int n, PLASMA_Complex64_t *A, int lda, int gM, int m0, int n0, unsigned long long int seed ),
          (m, n, A, lda, gM, m0, n0, seed) )
FUNCTION_TYPE( CORE_zpltmg, PLRNT, int ,
          ( PLASMA_enum mtxtype, int M, int N, PLASMA_Complex64_t *A, int LDA, int gM, int gN, int m0, int n0, unsigned long long int seed ),
          (mtxtype, M, N, A, LDA, gM, gN, m0, n0, seed) )
FUNCTION_TYPE( CORE_zpltmg_chebvand, PLRNT, int ,
          ( int M, int N, PLASMA_Complex64_t *A, int LDA, int gN, int m0, int n0, PLASMA_Complex64_t *W ),
          (M, N, A, LDA, gN, m0, n0, W) )
FUNCTION_TYPE( CORE_zpltmg_circul, PLRNT, int ,
          ( int M, int N, PLASMA_Complex64_t *A, int LDA, int gM, int m0, int n0, const PLASMA_Complex64_t *V ),
          (M, N, A, LDA, gM, m0, n0, V) )
FUNCTION_VOID( CORE_zpltmg_condexq, PLRNT, void ,
          ( int M, int N, PLASMA_Complex64_t *Q, int LDQ ),
          (M, N, Q, LDQ) )
FUNCTION_VOID( CORE_zpltmg_fiedler, PLRNT, void ,
          ( int M, int N, const PLASMA_Complex64_t *X, int incX, const PLASMA_Complex64_t *Y, int incY, PLASMA_Complex64_t *A, int LDA ),
          (M, N, X, incX, Y, incY, A, LDA) )
FUNCTION_TYPE( CORE_zpltmg_hankel, PLRNT, int ,
          ( PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t *A, int LDA, int m0, int n0, int nb, const PLASMA_Complex64_t *V1, const PLASMA_Complex64_t *V2 ),
          (uplo, M, N, A, LDA, m0, n0, nb, V1, V2) )
FUNCTION_VOID( CORE_zpltmg_toeppd1, PLRNT, void ,
          ( int gM, int m0, int M, PLASMA_Complex64_t *W, unsigned long long int seed ),
          (gM, m0, M, W, seed) )
FUNCTION_VOID( CORE_zpltmg_toeppd2, PLRNT, void ,
          ( int M, int N, int K, int m0, int n0, const PLASMA_Complex64_t *W, PLASMA_Complex64_t *A, int LDA ),
          (M, N, K, m0, n0, W, A, LDA) )
FUNCTION_VOID( CORE_zpotrf, POTRF, void ,
          (PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, int *INFO),
          (uplo, N, A, LDA, INFO) )
FUNCTION_VOID( CORE_zsetvar, SETVAR, void ,
          (const PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *x),
          (alpha, x) )
FUNCTION_VOID( CORE_zshiftw, SHIFTW, void ,
          (int s, int cl, int m, int n, int L, PLASMA_Complex64_t *A, PLASMA_Complex64_t *W) ,
          (s, cl, m, n, L, A, W)  )
FUNCTION_VOID( CORE_zshift, SHIFT, void ,
          (int s, int m, int n, int L, PLASMA_Complex64_t *A) ,
          (s, m, n, L, A)  )
FUNCTION_TYPE( CORE_zssssm, SSSSM, int ,
          (int M1, int N1, int M2, int N2, int K, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, const PLASMA_Complex64_t *L1, int LDL1, const PLASMA_Complex64_t *L2, int LDL2, const int *IPIV),
          (M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, L1, LDL1, L2, LDL2, IPIV) )
FUNCTION_VOID( CORE_zswpab, SWPAB, void ,
          (int i, int n1, int n2, PLASMA_Complex64_t *A, PLASMA_Complex64_t *work) ,
          (i, n1, n2, A, work)  )
FUNCTION_VOID( CORE_zsymm, SYMM, void ,
          (PLASMA_enum side, PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC),
          (side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC) )
FUNCTION_VOID( CORE_zsyr2k, SYR2K, void ,
          (PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *B, int LDB, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, B, LDB, beta, C, LDC) )
FUNCTION_VOID( CORE_zsyrk, SYRK, void ,
          (PLASMA_enum uplo, PLASMA_enum trans, int N, int K, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC),
          (uplo, trans, N, K, alpha, A, LDA, beta, C, LDC) )
FUNCTION_TYPE( CORE_zsyssq, LASSQ, int ,
          (PLASMA_enum uplo, int N, const PLASMA_Complex64_t *A, int LDA, double *scale, double *sumsq),
          (uplo, N, A, LDA, scale, sumsq) )
FUNCTION_VOID( CORE_ztrasm, ASUM, void ,
          (PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int M, int N, const PLASMA_Complex64_t *A, int lda, double *work),
          (storev, uplo, diag, M, N, A, lda, work) )
FUNCTION_VOID( CORE_ztrdalg1, TRDALG, void ,
          ( int n, int nb, PLASMA_Complex64_t *A, int lda, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU, int Vblksiz, int wantz, int i, int sweepid, int m, int grsiz, PLASMA_Complex64_t *work),
          (n, nb, A, lda, V, TAU, Vblksiz, wantz, i, sweepid, m, grsiz, work) )
FUNCTION_VOID( CORE_ztrmm, TRMM, void ,
          (PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int M, int N, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB),
          (side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB) )
FUNCTION_VOID( CORE_ztrsm, TRSM, void ,
          (PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int M, int N, PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *B, int LDB),
          (side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB) )
FUNCTION_TYPE( CORE_ztrssq, LASSQ, int ,
          (PLASMA_enum uplo, PLASMA_enum diag, int M, int N, const PLASMA_Complex64_t *A, int LDA, double *scale, double *sumsq),
          (uplo, diag, M, N, A, LDA, scale, sumsq) )
FUNCTION_VOID( CORE_ztrtri, TRTRI, void ,
          (PLASMA_enum uplo, PLASMA_enum diag, int N, PLASMA_Complex64_t *A, int LDA, int *info),
          (uplo, diag, N, A, LDA, info) )
FUNCTION_TYPE( CORE_ztslqt, TSLQT, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_ztsmlq, TSMLQ, int ,
          (PLASMA_enum side, PLASMA_enum trans, int M1, int N1, int M2, int N2, int K, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, const PLASMA_Complex64_t *V, int LDV, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_ztsmlq_corner, TSMLQ, int ,
          ( int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb, PLASMA_Complex64_t *A1, int lda1, PLASMA_Complex64_t *A2, int lda2, PLASMA_Complex64_t *A3, int lda3, const PLASMA_Complex64_t *V, int ldv, const PLASMA_Complex64_t *T, int ldt, PLASMA_Complex64_t *WORK, int ldwork),
          (m1, n1, m2, n2, m3, n3, k, ib, nb, A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_ztsmlq_hetra1, TSMLQ, int ,
          ( PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, PLASMA_Complex64_t *A1, int lda1, PLASMA_Complex64_t *A2, int lda2, const PLASMA_Complex64_t *V, int ldv, const PLASMA_Complex64_t *T, int ldt, PLASMA_Complex64_t *WORK, int ldwork),
          (side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_ztsmqr, TSMQR, int ,
          (PLASMA_enum side, PLASMA_enum trans, int M1, int N1, int M2, int N2, int K, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, const PLASMA_Complex64_t *V, int LDV, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_ztsmqr_corner, TSMQR, int ,
          ( int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb, PLASMA_Complex64_t *A1, int lda1, PLASMA_Complex64_t *A2, int lda2, PLASMA_Complex64_t *A3, int lda3, const PLASMA_Complex64_t *V, int ldv, const PLASMA_Complex64_t *T, int ldt, PLASMA_Complex64_t *WORK, int ldwork),
          (m1, n1, m2, n2, m3, n3, k, ib, nb, A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_ztsmqr_hetra1, TSMQR, int ,
          ( PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, PLASMA_Complex64_t *A1, int lda1, PLASMA_Complex64_t *A2, int lda2, const PLASMA_Complex64_t *V, int ldv, const PLASMA_Complex64_t *T, int ldt, PLASMA_Complex64_t *WORK, int ldwork),
          (side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork) )
FUNCTION_TYPE( CORE_ztsqrt, TSQRT, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_ztstrf, TSTRF, int ,
          (int M, int N, int IB, int NB, PLASMA_Complex64_t *U, int LDU, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *L, int LDL, int *IPIV, PLASMA_Complex64_t *WORK, int LDWORK, int *INFO),
          (M, N, IB, NB, U, LDU, A, LDA, L, LDL, IPIV, WORK, LDWORK, INFO) )
FUNCTION_TYPE( CORE_zttlqt, TTLQT, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_zttmlq, TTMLQ, int ,
          (PLASMA_enum side, PLASMA_enum trans, int M1, int N1, int M2, int N2, int K, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, const PLASMA_Complex64_t *V, int LDV, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_zttmqr, TTMQR, int ,
          (PLASMA_enum side, PLASMA_enum trans, int M1, int N1, int M2, int N2, int K, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, const PLASMA_Complex64_t *V, int LDV, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, M1, N1, M2, N2, K, IB, A1, LDA1, A2, LDA2, V, LDV, T, LDT, WORK, LDWORK) )
FUNCTION_TYPE( CORE_zttqrt, TTQRT, int ,
          (int M, int N, int IB, PLASMA_Complex64_t *A1, int LDA1, PLASMA_Complex64_t *A2, int LDA2, PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK),
          (M, N, IB, A1, LDA1, A2, LDA2, T, LDT, TAU, WORK) )
FUNCTION_TYPE( CORE_zunmlq, UNMLQ, int ,
          (PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, int IB, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *C, int LDC, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK) )
FUNCTION_TYPE( CORE_zunmqr, UNMQR, int ,
          (PLASMA_enum side, PLASMA_enum trans, int M, int N, int K, int IB, const PLASMA_Complex64_t *A, int LDA, const PLASMA_Complex64_t *T, int LDT, PLASMA_Complex64_t *C, int LDC, PLASMA_Complex64_t *WORK, int LDWORK),
          (side, trans, M, N, K, IB, A, LDA, T, LDT, C, LDC, WORK, LDWORK) )

