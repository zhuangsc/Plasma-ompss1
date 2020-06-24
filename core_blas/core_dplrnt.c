/**
 *
 * @file core_dplrnt.c
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
 * @generated d Tue Jan  7 11:44:47 2014
 *
 **/
#include "common.h"
#include "random.h"

#define REAL
#undef COMPLEX

#ifdef COMPLEX
#define NBELEM   2
#else
#define NBELEM   1
#endif

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dplrnt generates a random tile.
 *
 *******************************************************************************
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
#pragma weak CORE_dplrnt = PCORE_dplrnt
#define CORE_dplrnt PCORE_dplrnt
#endif
void CORE_dplrnt( int m, int n, double *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed )
{
    double *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)gM;

    for (j=0; j<n; ++j ) {
        ran = Rnd64_jump( NBELEM*jump, seed );
        for (i = 0; i < m; ++i) {
            *tmp = 0.5f - ran * RndF_Mul;
            ran  = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
            *tmp += I*(0.5f - ran * RndF_Mul);
            ran   = Rnd64_A * ran + Rnd64_C;
#endif
            tmp++;
        }
        tmp  += lda-i;
        jump += gM;
    }
}
