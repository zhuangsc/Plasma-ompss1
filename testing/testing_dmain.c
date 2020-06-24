/**
 *
 * @file testing_dmain.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:45:17 2014
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <plasma.h>
#undef COMPLEX
#define REAL
#undef SINGLE
#define DOUBLE
#include "testing_dmain.h"

int   IONE     = 1;
int   ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */

int   format[6] = { PlasmaRM, PlasmaCM, PlasmaCCRB, PlasmaCRRB, PlasmaRCRB, PlasmaRRRB };
int   side[2]   = { PlasmaLeft,    PlasmaRight };
int   uplo[2]   = { PlasmaUpper,   PlasmaLower };
int   diag[2]   = { PlasmaNonUnit, PlasmaUnit  };
int   trans[3]  = { PlasmaNoTrans, PlasmaTrans, PlasmaTrans };
int   itype[3]  = { 1, 2, 3 };
int   storev[2] = { PlasmaRowwise, PlasmaColumnwise };
int   norm[4]   = { PlasmaMaxNorm, PlasmaOneNorm, PlasmaInfNorm, PlasmaFrobeniusNorm };


char *formatstr[6]= { "RM", "CM", "CCRB", "CRRB", "RCRB", "RRRB"};
char *sidestr[2]  = { "Left ", "Right" };
char *uplostr[2]  = { "Upper", "Lower" };
char *diagstr[2]  = { "NonUnit", "Unit   " };
char *transstr[3] = { "N", "T", "H" };
char *itypestr[3] = { "inv(U')xAxinv(U) or inv(L)xAxinv(L')", "UxAxU' or L'xAxL", "UxAxU' or L'xAxL" };
char *storevstr[2]= { "Rowwise", "Columnwise" };
char *normstr[4]  = { "Max", "One", "Inf", "Fro" };

#define map_cm(m, n, i, j)   ((i) + (j) * (m))
#define map_rm(m, n, i, j)   ((i) * (n) + (j))

int map_CM(  int m, int n, int mb, int nb, int i, int j) { return map_cm(m, n, i, j); }
int map_RM(  int m, int n, int mb, int nb, int i, int j) { return map_rm(m, n, i, j); }
int map_CCRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_cm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_cm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n0 + ( (j/nb) * (nb*(m%mb)) )       + map_cm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m*n0  + ( (i/mb) * (mb*(n%nb)) )       + map_cm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_cm( m%mb, n%nb, i%mb, j%nb) );
}

int map_CRRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_cm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_rm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n0 + ( (j/nb) * (nb*(m%mb)) )       + map_rm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m*n0  + ( (i/mb) * (mb*(n%nb)) )       + map_rm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_rm( m%mb, n%nb, i%mb, j%nb) );
}

int map_RCRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_rm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_cm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n  + ( (j/nb) * (nb*(m%mb)) )       + map_cm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m0*n0 + ( (i/mb) * (mb*(n%nb)) )       + map_cm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_cm( m%mb, n%nb, i%mb, j%nb) );
}

int map_RRRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_rm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_rm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n  + ( (j/nb) * (nb*(m%mb)) )       + map_rm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m0*n0 + ( (i/mb) * (mb*(n%nb)) )       + map_rm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_rm( m%mb, n%nb, i%mb, j%nb) );
}

double gtime()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec * 1e3 + 1e-3 * tp.tv_usec;
}

void dump_info(char* lname, int nthreads, double elapsed, int rep)
{
	FILE *log;
	log = fopen(lname, "w");
	fprintf(log, "time: %.5f ms\n", elapsed);
	fprintf(log, "reps: %d\n", rep);
	fprintf(log, "threads: %d\n", nthreads);
	fclose(log);
}


void *formatmap[6] = {  map_RM, map_CM, map_CCRB, map_CRRB, map_RCRB, map_RRRB };

int main (int argc, char **argv)
{
    int ncores, sched;
    int info;
	int runtime;
    char func[32];

    /* Check for number of arguments*/
    if ( argc < 4) {
        printf(" Proper Usage is : ./dtesting ncores sched FUNC ...\n"
               "   - ncores : number of cores \n"
               "   - sched  : 0 for static, 1 for dynamic\n"
               "   - RUN    : 0 for QUARK, 1 for OmpSs\n"
               "   - FUNC   : name of function to test\n");
        exit(1);
    }

    sscanf( argv[1], "%d", &ncores );
    sscanf( argv[2], "%d", &sched  );
    sscanf( argv[3], "%d", &runtime );
    sscanf( argv[4], "%s", func   );

    PLASMA_Runtime_Init(ncores, runtime);

    if ( sched == 0 )
        PLASMA_Set( PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );
    else
        PLASMA_Set( PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );

	if ( runtime == 0 )
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
	else 
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);

    argc -= 5;
    argv += 5;
    info  = 0;

    /*
     * Norms
     */
    if ( strcmp(func, "LANGE") == 0 ) {
        info = testing_dlange( argc, argv );
    /*
     * Blas Level 3
     */
    } else if ( strcmp(func, "GEMM") == 0 ) {
        info = testing_dgemm( argc, argv );
#ifdef COMPLEX
    } else if ( strcmp(func, "HEMM") == 0 ) {
        info = testing_dsymm( argc, argv );
    } else if ( strcmp(func, "HERK") == 0 ) {
        info = testing_dsyrk( argc, argv );
    } else if ( strcmp(func, "HER2K") == 0 ) {
      info = testing_dsyr2k( argc, argv );
#endif
    } else if ( strcmp(func, "SYMM") == 0 ) {
        info = testing_dsymm( argc, argv );
    } else if ( strcmp(func, "SYRK") == 0 ) {
        info = testing_dsyrk( argc, argv );
        } else if ( strcmp(func, "SYR2K") == 0 ) {
                info = testing_dsyr2k( argc, argv );
    } else if ( strcmp(func, "TRMM") == 0 ) {
        info = testing_dtrmm( argc, argv );
    } else if ( strcmp(func, "TRSM") == 0 ) {
        info = testing_dtrsm( argc, argv );
    } else if ( strcmp(func, "PEMV") == 0 ) {
        info = testing_dpemv( argc, argv );
    /*
     * Linear system
     */
    } else if ( strcmp(func, "POSV") == 0 ) {
        info = testing_dposv( argc, argv );
    } else if ( strcmp(func, "GELS") == 0 ) {
        info = testing_dgels( argc, argv );
    } else if ( strcmp(func, "GESV_INCPIV") == 0 ) {
        info = testing_dgesv_incpiv( argc, argv );
    } else if ( strcmp(func, "GESV") == 0 ) {
        info = testing_dgesv( argc, argv );
    /*
     * Matrix inversion
     */
    } else if ( strcmp(func, "POTRI") == 0 ) {
        info = testing_dpotri( argc, argv );
    } else if ( strcmp(func, "GETRI") == 0 ) {
        info = testing_dgetri( argc, argv );
    /*
     * Eigenvalue Problems
     */
    } else if ( strcmp(func, "HETRD") == 0 ) {
      info = testing_dsytrd( argc, argv );
    } else if ( strcmp(func, "HEEV") == 0 ) {
      info = testing_dsyev( argc, argv );
    } else if ( strcmp(func, "HEEVD") == 0 ) {
      info = testing_dsyevd( argc, argv );
    } else if ( strcmp(func, "HEEVR") == 0 ) {
      info = testing_dsyevr( argc, argv );
    } else if ( strcmp(func, "HEGV") == 0 ) {
      info = testing_dsygv( argc, argv );
    } else if ( strcmp(func, "HEGVD") == 0 ) {
      info = testing_dsygv( argc, argv );
    } else if ( strcmp(func, "HEGST") == 0 ) {
      info = testing_dsygst( argc, argv );
    /*
     * Singular Value Decomposition
     */
    } else if ( strcmp(func, "GEBRD") == 0 ) {
      info = testing_dgebrd( argc, argv );
    } else if ( strcmp(func, "GESVD") == 0 ) {
      info = testing_dgesvd( argc, argv );
    } else if ( strcmp(func, "GESDD") == 0 ) {
      info = testing_dgesdd( argc, argv );
#ifdef DOUBLE
    /*
     * Mixed precision
     */
    } else if ( strcmp(func, "SPOSV") == 0 ) {
        info = testing_dsposv( argc, argv );
    } else if ( strcmp(func, "SGESV") == 0 ) {
        info = testing_dsgesv( argc, argv );
    } else if ( strcmp(func, "SUNGESV") == 0 ) {
        info = testing_dsungesv( argc, argv );
#endif
    /* Layout Transformation */
    } else if ( strcmp(func, "GECFI") == 0 ) {
        info = testing_dgecfi( argc, argv );
    } else if ( strcmp(func, "GETMI") == 0 ) {
        info = testing_dgetmi( argc, argv );
	} else if ( strcmp(func, "POTRF") == 0 ) {
		info = testing_dpotrf( argc, argv );
	} else if ( strcmp(func, "GEQRF") == 0 ) {
		info = testing_dgeqrf(argc, argv );
	} else if ( strcmp(func, "GETRF") == 0 ) {
		info = testing_dgetrf(argc, argv);
	} else if ( strcmp(func, "SPDSOLV") == 0 ) {
		info = testing_dspdsolv(argc, argv);
	} else if ( strcmp(func, "GEMM_TIMING") == 0 ) {
		info = timing_dgemm(argc, argv);
	} else if ( strcmp(func, "POTRF_TIMING") == 0 ) {
		info = timing_dpotrf(argc, argv);
	} else if ( strcmp(func, "SYRK_TIMING") == 0 ) {
		info = timing_dsyrk(argc, argv);
    } else {
        fprintf(stderr, "Function unknown\n");
    }

    if ( info == -1 ) {
        printf( "TESTING %s FAILED : incorrect number of arguments\n", func);
    } else if ( info == -2 ) {
        printf( "TESTING %s FAILED : not enough memory\n", func);
    }

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}
