/**
 *
 * @generated c Tue Jan  7 11:45:26 2014
 *
 **/
#ifndef CAUXILIARY_H
#define CAUXILIARY_H

int    c_check_orthogonality   (int M, int N, int LDQ, PLASMA_Complex32_t *Q);
int    c_check_QRfactorization (int M, int N, PLASMA_Complex32_t *A1, PLASMA_Complex32_t *A2, int LDA, PLASMA_Complex32_t *Q);
int    c_check_LLTfactorization(int N, PLASMA_Complex32_t *A1, PLASMA_Complex32_t *A2, int LDA, int uplo);
float c_check_gemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                   PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA,
                   PLASMA_Complex32_t *B, int LDB,
                   PLASMA_Complex32_t beta, PLASMA_Complex32_t *Cplasma,
                   PLASMA_Complex32_t *Cref, int LDC,
                   float *Cinitnorm, float *Cplasmanorm, float *Clapacknorm );

float c_check_trsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
           int M, int NRHS, PLASMA_Complex32_t alpha,
           PLASMA_Complex32_t *A, int LDA,
           PLASMA_Complex32_t *Bplasma, PLASMA_Complex32_t *Bref, int LDB,
           float *Binitnorm, float *Bplasmanorm, float *Blapacknorm );

float c_check_solution(int M, int N, int NRHS,
                      PLASMA_Complex32_t *A1, int LDA,
                      PLASMA_Complex32_t *B1, PLASMA_Complex32_t *B2, int LDB,
                      float *anorm, float *bnorm, float *xnorm);

#endif /* CAUXILIARY_H */
