/**
 *
 * @generated s Tue Jan  7 11:45:26 2014
 *
 **/
#ifndef SAUXILIARY_H
#define SAUXILIARY_H

int    s_check_orthogonality   (int M, int N, int LDQ, float *Q);
int    s_check_QRfactorization (int M, int N, float *A1, float *A2, int LDA, float *Q);
int    s_check_LLTfactorization(int N, float *A1, float *A2, int LDA, int uplo);
float s_check_gemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                   float alpha, float *A, int LDA,
                   float *B, int LDB,
                   float beta, float *Cplasma,
                   float *Cref, int LDC,
                   float *Cinitnorm, float *Cplasmanorm, float *Clapacknorm );

float s_check_trsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
           int M, int NRHS, float alpha,
           float *A, int LDA,
           float *Bplasma, float *Bref, int LDB,
           float *Binitnorm, float *Bplasmanorm, float *Blapacknorm );

float s_check_solution(int M, int N, int NRHS,
                      float *A1, int LDA,
                      float *B1, float *B2, int LDB,
                      float *anorm, float *bnorm, float *xnorm);

#endif /* SAUXILIARY_H */
