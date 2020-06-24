/**
 *
 * @file core_blas-gpu.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:39 2014
 *
 **/
#ifndef _PLASMA_CORE_BLASGPU_H_
#define _PLASMA_CORE_BLASGPU_H_

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of GPU kernels - alphabetical order
 **/

int cublasDpamm(cublasHandle_t handle, int op, cublasSideMode_t side, PLASMA_enum storev,
           int M, int N, int K, int L,
           const double *A1, int LDA1,
                 double *A2, int LDA2,
           const double *V, int LDV,
                 double *W, int LDW);

int cublasDparfb(cublasHandle_t handle, cublasSideMode_t side, cublasOperation_t trans, PLASMA_enum direct, PLASMA_enum storev,
            int M1, int N1, int M2, int N2, int K, int L,
                  double *A1, int LDA1,
                  double *A2, int LDA2,
            const double *V, int LDV,
            const double *T, int LDT,
                  double *WORK, int LDWORK);

int cublasDtsmqr(cublasHandle_t handle, cublasSideMode_t side, cublasOperation_t trans,
                int M1, int N1, int M2, int N2, int K, int IB,
                double *A1, int LDA1,
                double *A2, int LDA2,
                const double *V, int LDV,
                const double *T, int LDT,
                double *WORK, int LDWORK);

int cublasDlacpy(cudaStream_t stream, PLASMA_enum uplo, int m, int n,
                double *dA, int ldda,
                double *dB, int lddb);

#ifdef __cplusplus
}
#endif

#endif //_PLASMA_CORE_BLASGPU_H_
