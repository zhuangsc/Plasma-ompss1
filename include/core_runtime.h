#ifndef __RT_CORE_DGEMM__
#define __RT_CORE_DGEMM__

#include "quark.h"
#include "plasma.h"

void RT_CORE_dpltmg( Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum mtxtype, int m, int n, double *A, int lda,
                         int gM, int gN, int m0, int n0, unsigned long long int seed );

void RT_CORE_dplgsy( Quark *quark, Quark_Task_Flags *task_flags,
                        double bump, int m, int n, double *A, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed );

void RT_CORE_dlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda);

void RT_CORE_dtrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag,
                       int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);

void RT_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       const double *A, int lda,
                       double *B, int ldb);

void RT_CORE_dgeam(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, int transB,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *B, int ldb);

void RT_CORE_dscal(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, double *alpha, double *A, int lda);

void RT_CORE_dgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, int transB,
                      int m, int n, int k, int nb,
                      double alpha, double *A, int lda,
					  double *B, int ldb,
                      double beta, double *C, int ldc);

void RT_CORE_dgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, int transB,
                        int m, int n, int k, int nb,
                        double alpha, const double *A, int lda,
                        const double *B, int ldb,
                        double beta, double *C, int ldc);

void RT_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *C, int ldc);


void RT_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       double alpha, double *A, int lda,
                       double *B, int ldb,
                       double beta, double *C, int ldc);


void RT_CORE_dtrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb);

void RT_CORE_dtrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb);

void RT_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc);

void RT_CORE_dpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);

void RT_CORE_dlaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       double alpha, double beta,
                       double *A, int LDA);

void RT_CORE_dgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt);

void RT_CORE_dormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const double *A, int lda,
                       const double *T, int ldt,
                       double *C, int ldc);

void RT_CORE_dtsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);

void RT_CORE_dtsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt);

void RT_CORE_dlange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result);

void RT_CORE_dlange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, int M, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);

void RT_CORE_free(Quark *quark, Quark_Task_Flags *task_flags, void *A, int szeA);

void RT_CORE_dasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       const double *A, int lda, int szeA,
                       double *work, int szeW);

void RT_CORE_dasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          const double *A, int lda, int szeA,
                          double *work, int szeW, double *fake, int szeF);

void RT_CORE_dgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );

void RT_CORE_dplssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const double *A, double *result );

void CORE_dplssq(int m, double *A, double *result);

void RT_CORE_dgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);


void RT_LAPACKE_dgetrf_work(int m, int n, double *A, int lda, int *IPIV);

void RT_CORE_dgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);

void RT_CORE_dlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, double *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc);

void RT_CORE_dlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, double *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          double *fake1, int szefake1, int flag1,
                          double *fake2, int szefake2, int flag2);

void RT_CORE_dgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, int transB,
                         int m, int n, int k, int nb,
                         double alpha, const double *A, int lda,
                                                   const double *B, int ldb,
                         double beta, double *C, int ldc,
                         double *fake1, int szefake1, int flag1,
                         double *fake2, int szefake2, int flag2);

void RT_CORE_dsyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );

void RT_CORE_dgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, double *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);

void RT_CORE_dswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, double *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const double *Akk, int ldak);

void RT_CORE_dlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, double *Aij,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 double *fake1, int szefake1, int flag1,
                                 double *fake2, int szefake2, int flag2);

void RT_CORE_foo(double *A, int sze);

void RT_CORE_foo2(double *A, int szeA, double *B, int szeB);

#endif
