#include "runtime.h"

void RT_CORE_dscal(Quark *quark, Quark_Task_Flags *task_flags,
                      int m,  double *alpha, double *A, int lda)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
                cblas_dscal(m, 1.0/alpha[0], A, lda);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                //double *Adel = A + (m-256)*m;
                //#pragma omp task in([1]alpha) inout([m]Adel) no_copy_deps label(dscal)  
                #pragma omp target device (smp) copy_deps
                #pragma omp task in([1]alpha) inout([m]A) label(dscal)  
                cblas_dscal(m, 1.0/alpha[0], A, lda);
	}
}

void RT_CORE_dconv(Quark *quark, Quark_Task_Flags *task_flags,  
                      double *alpha, double *beta, double *result, double *e0)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
        #pragma omp target device (smp) copy_deps
        #pragma omp task in([1]alpha, [1]beta) inout([1]result) out([1]e0) label(scalars_for_normest) 
        {
                e0[0] = result[0];
                result[0] = alpha[0]/beta[0];
        }
	
}
void RT_CORE_conv( Quark *quark, Quark_Task_Flags *task_flags, double *e, double *e0, double *conv)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
        #pragma omp target device (smp) copy_deps
        #pragma omp task in([1]e, [1]e0) out([1]conv) label(conv_for_normest)  
        {
                conv[0] = abs( e[0] - e0[0]);
        }
	
}
void RT_CORE_tolconv( Quark *quark, Quark_Task_Flags *task_flags, double *tol, double *result, double *tolconv, double *conv)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
        #pragma omp target device (smp) copy_deps
        #pragma omp task in([1]tol, [1]result) out([1]tolconv, [1]conv) label(tolconv_for_normest)  
        {
                tolconv[0] = tol[0] * result[0];
                conv[0] = result[0];
        }
	
}

void RT_CORE_dgemm_2norm( Quark *quark, Quark_Task_Flags *task_flags, PLASMA_enum transA, 
                          PLASMA_enum transB, double alpha, PLASMA_desc A, double *B, 
                          int szeB, double beta, double *C, int szeC)
{
     plasma_context_t *plasma;
     plasma = plasma_context_self();
     PLASMA_desc *descB2;  PLASMA_Desc_Create( &descB2 , B, PlasmaRealDouble, A.nb, 1, A.nb, A.n, 1, 0, 0, A.n, 1);
     PLASMA_desc descB = plasma_desc_submatrix(*descB2, 0, 0, A.n, 1);
     PLASMA_desc *descC2; PLASMA_Desc_Create( &descC2, C, PlasmaRealDouble, A.nb, 1, A.nb, A.n, 1, 0, 0, A.n, 1);
     PLASMA_desc descC = plasma_desc_submatrix(*descC2, 0, 0, A.n, 1);
     int nnz, i;
     int IONE=1;
     int ISEED[4] = {0,0,0,1};
      
     //double *Cdel = BLKADDR(descC, double,descC.mt-1, descC.nt-1); 
     //double *Bdel = BLKADDR(descB, double, descC.mt-1, descC.nt-1); 
     //#pragma omp task in([descC.nb*descC.nb]Bdel) inout([szeC]C) no_copy_deps label(dgemm_2norm) 
     #pragma omp target device (smp) copy_deps
     #pragma omp task in([szeB]B) inout([szeC]C) label(dgemm_2norm) 
     RT_CORE_dgemm_2norm1( transA, transB, alpha, A, B, descB, 
                          szeB, beta, C, descC, szeC);
     #pragma omp target device (smp) copy_deps
     #pragma omp task inout([szeC]C) label(dgemm_2norm) 
     {
        for(i= 0 ; i < szeC ; i++ ){
           if ( C[i] != 0)
               nnz = nnz + 1;
        }
        if ( nnz == 0){
            LAPACKE_dlarnv_work(IONE, ISEED, szeC, C);
        }
     }

     //#pragma omp task in([descC.nb*descC.nb]Cdel) inout([szeB]B) no_copy_deps label(dgemm_2norm) 
     #pragma omp target device (smp) copy_deps
     #pragma omp task in([szeC]C) inout([szeB]B) label(dgemm_2norm)
     RT_CORE_dgemm_2norm2( transA, transB, alpha, A, C, descC, 
                          szeC, beta, B, descB, szeB);

}



#define Ad(m, n) BLKADDR(A, double, m, n)
#define descB(m, n) BLKADDR(descB, double, m, n)
#define descC(m, n) BLKADDR(descC, double, m, n)
void RT_CORE_dgemm_2norm1( PLASMA_enum transA, PLASMA_enum transB, double alpha, PLASMA_desc A, 
                           double *B, PLASMA_desc descB,
                          int szeB, double beta, 
                          double *C, PLASMA_desc descC, int szeC)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
        int ldam, ldak, ldbn, ldbk, ldcm, m, n, k;
        int tempmm, tempnn, tempkn, tempkm;

        double zbeta;
        double zone = (double)1.0;
           for (m = 0; m < descC.mt; m++) {
              tempmm = m == descC.mt-1 ? descC.m-m*descC.mb : descC.mb;
              ldcm = BLKLDD(descC, m);
               for (n = 0; n < descC.nt; n++) {
                    tempnn = n == descC.nt-1 ? descC.n-n*descC.nb : descC.nb;
                    ldam = BLKLDD(A, m);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        ldbk = BLKLDD(descB, k);
                        zbeta = k == 0 ? beta : zone;
                        CORE_dgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkn, 
                            1.0, Ad(m, k), ldam,  
                                 descB(k, n), ldbk,  
                            zbeta, descC(m, n), ldcm); 
                    }
               }
           }
}
void RT_CORE_dgemm_2norm2( PLASMA_enum transA, PLASMA_enum transB, double alpha, PLASMA_desc A, 
                           double *C, PLASMA_desc descC,
                          int szeC, double beta, 
                          double *B, PLASMA_desc descB, int szeB)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
        int ldam, ldak, ldbn, ldbk, ldcm, m, n, k;
        int tempmm, tempnn, tempkn, tempkm;

        double zbeta;
        double zone = (double)1.0;
           for (m = 0; m < descC.mt; m++) {
              tempmm = m == descC.mt-1 ? descC.m-m*descC.mb : descC.mb;
              ldcm = BLKLDD(descC, m);
               for (n = 0; n < descC.nt; n++) {
                    tempnn = n == descC.nt-1 ? descC.n-n*descC.nb : descC.nb;
                    for (k = 0; k < A.nt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(descC, k);
                        zbeta = k == 0 ? beta : zone;
                        CORE_dgemm(
                            PlasmaTrans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm, 
                            1.0, Ad(k, m), ldak,  
                                 descC(k, n), ldbk,  
                            zbeta, descB(m, n), ldcm); 
                    }
               }
           }
}
