/**
 *
 * @generated d Tue Jan  7 11:45:26 2014
 *
 **/
#ifndef DAUXILIARY_H
#define DAUXILIARY_H

int    d_check_orthogonality   (int M, int N, int LDQ, double *Q);
int    d_check_QRfactorization (int M, int N, double *A1, double *A2, int LDA, double *Q);
int    d_check_LLTfactorization(int N, double *A1, double *A2, int LDA, int uplo);
double d_check_gemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                   double alpha, double *A, int LDA,
                   double *B, int LDB,
                   double beta, double *Cplasma,
                   double *Cref, int LDC,
                   double *Cinitnorm, double *Cplasmanorm, double *Clapacknorm );

double d_check_trsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
           int M, int NRHS, double alpha,
           double *A, int LDA,
           double *Bplasma, double *Bref, int LDB,
           double *Binitnorm, double *Bplasmanorm, double *Blapacknorm );

double d_check_solution(int M, int N, int NRHS,
                      double *A1, int LDA,
                      double *B1, double *B2, int LDB,
                      double *anorm, double *bnorm, double *xnorm);

#endif /* DAUXILIARY_H */
