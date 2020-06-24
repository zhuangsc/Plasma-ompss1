/**
 *
 * @file core_zplghe.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 **/
#include "common.h"
#include "random.h"

#define COMPLEX
#undef REAL

#ifdef COMPLEX
#define NBELEM   2
#else
#define NBELEM   1
#endif

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zplgsy generates an hermitian matrix.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *         Scalar added to the diagonal of the full Matrix A to make it diagonal
 *         dominant.
 *
 * @param[in] m
 *         The number of rows of the tile A. m >= 0.
 *
 * @param[in] n
 *         The number of columns of the tile A. n >= 0.
 *
 * @param[in,out] A
 *         On entry, the m-by-n tile to be initialized.
 *         On exit, the tile initialized in the mtxtype format.
 *
 * @param[in] lda
 *         The leading dimension of the tile A. lda >= max(1,m).
 *
 * @param[in] gM
 *         The global number of rows of the full matrix, A is belonging to. gM >= (m0+M).
 *
 * @param[in] m0
 *         The index of the first row of tile A in the full matrix. m0 >= 0.
 *
 * @param[in] n0
 *         The index of the first column of tile A in the full matrix. n0 >= 0.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zplghe = PCORE_zplghe
#define CORE_zplghe PCORE_zplghe
#endif
void CORE_zplghe( double bump, int m, int n, PLASMA_Complex64_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed )
{
    PLASMA_Complex64_t *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)gM;

    /*
     * Tile diagonal
     */
    if ( m0 == n0 ) {
        for (j = 0; j < n; j++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (i = j; i < m; i++) {
                *tmp = 0.5f - ran * RndF_Mul;
                ran  = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
                *tmp += I*(0.5f - ran * RndF_Mul);
                ran   = Rnd64_A * ran + Rnd64_C;
#endif
                tmp++;
            }
            tmp  += (lda - i + j + 1);
            jump += gM + 1;
        }

        for (j = 0; j < n; j++) {
#ifdef COMPLEX
            A[j+j*lda] += bump - I*cimag( A[j+j*lda] );
#else
            A[j+j*lda] += bump;
#endif

            for (i=0; i<j; i++) {
                A[lda*j+i] = conj( A[lda*i+j] );
            }
        }
    }
    /*
     * Lower part
     */
    else if ( m0 > n0 ) {
        for (j = 0; j < n; j++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (i = 0; i < m; i++) {
                *tmp = 0.5f - ran * RndF_Mul;
                ran  = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
                *tmp += I*(0.5f - ran * RndF_Mul);
                ran   = Rnd64_A * ran + Rnd64_C;
#endif
                tmp++;
            }
            tmp  += (lda - i);
            jump += gM;
        }
    }
    /*
     * Upper part
     */
    else if ( m0 < n0 ) {
        /* Overwrite jump */
        jump = n0 + m0 * gM;

        for (i = 0; i < m; i++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (j = 0; j < n; j++) {
                A[j*lda+i] = 0.5f - ran * RndF_Mul;
                ran = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
                A[j*lda+i] -= I*(0.5f - ran * RndF_Mul);
                ran = Rnd64_A * ran + Rnd64_C;
#endif
            }
            jump += gM;
        }
    }
}
