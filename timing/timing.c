/**
 *
 * @file timing.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author ???
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 *
 **/

#if defined( _WIN32 ) || defined( _WIN64 )
#define int64_t __int64
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef PLASMA_EZTRACE
#include <eztrace.h>
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#include <time.h>
#include <sys/timeb.h>
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    FILETIME         ft;
    unsigned __int64 tmpres = 0;
    static int       tzflag;

    if (NULL != tv)
        {
            GetSystemTimeAsFileTime(&ft);
            tmpres |=  ft.dwHighDateTime;
            tmpres <<= 32;
            tmpres |=  ft.dwLowDateTime;

            /*converting file time to unix epoch*/
            tmpres /= 10;  /*convert into microseconds*/
            tmpres -= DELTA_EPOCH_IN_MICROSECS;

            tv->tv_sec  = (long)(tmpres / 1000000UL);
            tv->tv_usec = (long)(tmpres % 1000000UL);
        }
    if (NULL != tz)
        {
            if (!tzflag)
                {
                    _tzset();
                    tzflag++;
                }
            tz->tz_minuteswest = _timezone / 60;
            tz->tz_dsttime     = _daylight;
        }
    return 0;
}

#else  /* Non-Windows */
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <cblas.h>
#include <lapacke.h>
#include <plasma.h>
#include <core_blas.h>
#include "flops.h"
#include "timing.h"
#include "auxiliary.h"
#include <omp.h>

static int RunTest(int *iparam, _PREC *dparam, double *t_);

double cWtime(void);

int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

/*
 * struct timeval {time_t tv_sec; suseconds_t tv_usec;};
 */
double cWtime(void)
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

static int
Test(int64_t n, int *iparam) {
    int      i, j, iter;
    int      thrdnbr, niter;
    int64_t  M, N, K, NRHS;
    double  *t;
    _PREC    eps = _LAMCH( 'e' );
    _PREC    dparam[IPARAM_DNBPARAM];
    double   fmuls, fadds, fp_per_mul, fp_per_add;
    double   sumgf, sumgf2, sumt, sd, flops, gflops;
    char    *s;
    char    *env[] = {
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "ACML_NUM_THREADS",
        "ATLAS_NUM_THREADS",
        "BLAS_NUM_THREADS", ""
    };
    int gnuplot = 0;

    memset( &dparam, 0, IPARAM_DNBPARAM * sizeof(_PREC) );

    thrdnbr = iparam[IPARAM_THRDNBR];
    niter   = iparam[IPARAM_NITER];

    M    = iparam[IPARAM_M];
    N    = iparam[IPARAM_N];
    K    = iparam[IPARAM_K];
    NRHS = K;
    (void)M;(void)N;(void)K;(void)NRHS;

    if ( (n < 0) || (thrdnbr < 0 ) ) {
        if (gnuplot) {
            printf( "set title '%d_NUM_THREADS: ", thrdnbr );
            for (i = 0; env[i][0]; ++i) {
                s = getenv( env[i] );

                if (i) printf( " " ); /* separating space */

                for (j = 0; j < 5 && env[i][j] && env[i][j] != '_'; ++j)
                    printf( "%c", env[i][j] );

                if (s)
                    printf( "=%s", s );
                else
                    printf( "->%s", "?" );
            }
            printf( "'\n" );
            printf( "%s\n%s\n%s\n%s\n%s%s%s\n",
                    "set xlabel 'Matrix size'",
                    "set ylabel 'Gflop/s'",
                    "set key bottom",
                    gnuplot > 1 ? "set terminal png giant\nset output 'timeplot.png'" : "",
                    "plot '-' using 1:5 title '", _NAME, "' with linespoints" );
        }

        return 0;
    }

    printf( "%7d %7d %7d ", iparam[IPARAM_M], iparam[IPARAM_N], iparam[IPARAM_K] );
    fflush( stdout );

    t = (double*)malloc(niter*sizeof(double));
    memset(t, 0, niter*sizeof(double));

    if (sizeof(_TYPE) == sizeof(_PREC)) {
        fp_per_mul = 1;
        fp_per_add = 1;
    } else {
        fp_per_mul = 6;
        fp_per_add = 2;
    }

    fadds = (double)(_FADDS);
    fmuls = (double)(_FMULS);
    flops = 1e-9 * (fmuls * fp_per_mul + fadds * fp_per_add);
    gflops = 0.0;

    if ( iparam[IPARAM_WARMUP] ) {
        RunTest( iparam, dparam, &(t[0]));
    }

    sumgf  = 0.0;
    sumgf2 = 0.0;
    sumt   = 0.0;

    for (iter = 0; iter < niter; iter++)
    {
        if( iter == 0 ) {
          if ( iparam[IPARAM_TRACE] )
            iparam[IPARAM_TRACE] = 2;
          if ( iparam[IPARAM_DAG] )
            iparam[IPARAM_DAG] = 2;

          RunTest( iparam, dparam, &(t[iter]));

          iparam[IPARAM_TRACE] = 0;
          iparam[IPARAM_DAG] = 0;
        }
        else
            RunTest( iparam, dparam, &(t[iter]));

        gflops = flops / t[iter];

        sumt   += t[iter];
        sumgf  += gflops;
        sumgf2 += gflops*gflops;
    }

    gflops = sumgf / niter;
    sd = sqrt((sumgf2 - (sumgf*sumgf)/niter)/niter);

    if ( iparam[IPARAM_CHECK] )
        printf( "%9.3f %9.2f %9.2f %8.5e %8.5e %8.5e %8.5e %8.5e\n",
                sumt/niter, gflops, sd, dparam[IPARAM_RES], dparam[IPARAM_ANORM], dparam[IPARAM_XNORM], dparam[IPARAM_BNORM],
                dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ));
    else{
        printf( "%9.3f %9.2f %9.2f\n", sumt/niter, gflops, sd );
    }
    fflush( stdout );
    free(t);

    return 0;
}

static int
startswith(const char *s, const char *prefix) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}

static int
get_range(char *range, int *start_p, int *stop_p, int *step_p) {
    char *s, *s1, buf[21];
    int colon_count, copy_len, nbuf=20, n;
    int start=1000, stop=10000, step=1000;

    colon_count = 0;
    for (s = strchr( range, ':'); s; s = strchr( s+1, ':'))
        colon_count++;

    if (colon_count == 0) { /* No colon in range. */
        if (sscanf( range, "%d", &start ) < 1 || start < 1)
            return -1;
        step = start / 10;
        if (step < 1) step = 1;
        stop = start + 10 * step;

    } else if (colon_count == 1) { /* One colon in range.*/
        /* First, get the second number (after colon): the stop value. */
        s = strchr( range, ':' );
        if (sscanf( s+1, "%d", &stop ) < 1 || stop < 1)
            return -1;

        /* Next, get the first number (before colon): the start value. */
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start > stop || start < 1)
            return -1;

        /* Let's have 10 steps or less. */
        step = (stop - start) / 10;
        if (step < 1)
            step = 1;
    } else if (colon_count == 2) { /* Two colons in range. */
        /* First, get the first number (before the first colon): the start value. */
        s = strchr( range, ':' );
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start < 1)
            return -1;

        /* Next, get the second number (after the first colon): the stop value. */
        s1 = strchr( s+1, ':' );
        n = s1 - (s + 1);
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, s+1, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &stop ) < 1 || stop < start)
            return -1;

        /* Finally, get the third number (after the second colon): the step value. */
        if (sscanf( s1+1, "%d", &step ) < 1 || step < 1)
            return -1;
    } else

        return -1;

    *start_p = start;
    *stop_p = stop;
    *step_p = step;

    return 0;
}

static void
show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "\n"
            "  --[a]sync        Enable/Disable synchronous calls in wrapper function such as POTRI. (default: async)\n"
            "  --[no]check      Check result (default: nocheck)\n"
            "  --[no]inv        Check on inverse (default: noinv)\n"
            "  --[no]warmup     Perform a warmup run to pre-load libraries (default: warmup)\n"
#if defined(PLASMA_EZTRACE)
            "  --[no]trace      Enable/Disable trace generation (default: notrace)\n"
#endif
            "  --[no]dag        Enable/Disable DAG generation (default: nodag)\n"
            "                   Generates a dot_dag_file.dot.\n"
            "  --[no]atun       Activate autotuning (default: noatun) [Deprecated]\n"
            "  --inplace        Enable layout conversion inplace for lapack interface timers (default: enable)\n"
            "  --outplace       Enable layout conversion out of place for lapack interface timers (default: disable)\n"
            "\n"
            "  --n_range=R      Range of N values\n"
            "                   with R=Start:Stop:Step (default: 500:5000:500)\n"
            "  --m=X            dimension (M) of the matrices (default: N)\n"
            "  --k=X            dimension (K) of the matrices (default: 1)\n"
            "  --nrhs=X         Number of right-hand size (default: 1)\n"
            "  --nb=N           Nb size. (default: 128)\n"
            "  --ib=N           IB size. (default: 32)\n"
            "\n"
            "  --niter=N        Number of iterations performed for each test (default: 1)\n"
            "\n"
            "  --rhblk=N        If N > 0, enable Householder mode for QR and LQ factorization\n"
            "                   N is the size of each subdomain (default: 0)\n"
            "  --tntmode=x      Factorization used to select the pivots in the tournament pivoting algorithm (default: 1)\n"
            "                     1 - LU factorization with partial pivoting\n"
            "                     2 - QR rank revealing factorization\n"
            "  --tntsize=N      If N > 0, number of tiles grouped together during the tournament (default=4)\n"
            "\n"
            " Options specific to the conversion format timings xgetri and xgecfi:\n"
            "  --ifmt           Input format. (default: 0)\n"
            "  --ofmt           Output format. (default: 1)\n"
            "                   The possible values are:\n"
            "                     0 - PlasmaCM, Column major\n"
            "                     1 - PlasmaCCRB, Column-Colum rectangular block\n"
            "                     2 - PlasmaCRRB, Column-Row rectangular block\n"
            "                     3 - PlasmaRCRB, Row-Colum rectangular block\n"
            "                     4 - PlasmaRRRB, Row-Row rectangular block\n"
            "                     5 - PlasmaRM, Row Major\n"
            "  --thrdbypb       Number of threads per subproblem for inplace transformation (default: 1)\n"
            "\n");
}


static void
print_header(char *prog_name, int * iparam) {
    _PREC    eps = _LAMCH( 'e' );
    printf( "#\n"
            "# PLASMA %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            PLASMA_VERSION_MAJOR,
            PLASMA_VERSION_MINOR,
            PLASMA_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            iparam[IPARAM_NB],
            iparam[IPARAM_IB],
            eps );
    if  (iparam[IPARAM_CHECK] )
        printf( "#     M       N  K/NRHS   seconds   Gflop/s Deviation"
                "    ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)\n" );
    else
        printf( "#     M       N  K/NRHS   seconds   Gflop/s Deviation\n" );

    return;
}

static void
get_thread_count(int *thrdnbr) {
#if defined WIN32 || defined WIN64
    sscanf( getenv( "NUMBER_OF_PROCESSORS" ), "%d", thrdnbr );
#else
    *thrdnbr = sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

int
main(int argc, char *argv[]) {
    int i, m, mx, nx;
    int start =  500;
    int stop  = 5000;
    int step  =  500;
    int iparam[IPARAM_SIZEOF];

    memset(iparam, 0, IPARAM_SIZEOF*sizeof(int));

    iparam[IPARAM_THRDNBR       ] = 1;
    iparam[IPARAM_THRDNBR_SUBGRP] = 1;
    iparam[IPARAM_SCHEDULER     ] = 1;
    iparam[IPARAM_RUNTIME       ] = 2;
    iparam[IPARAM_M             ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_K             ] = 500;
    iparam[IPARAM_LDA           ] = 500;
    iparam[IPARAM_LDB           ] = 500;
    iparam[IPARAM_LDC           ] = 500;
    iparam[IPARAM_MB            ] = 128;
    iparam[IPARAM_NB            ] = 128;
    iparam[IPARAM_IB            ] = 32;
    iparam[IPARAM_NITER         ] = 1;
    iparam[IPARAM_WARMUP        ] = 1;
    iparam[IPARAM_CHECK         ] = 1;
    iparam[IPARAM_VERBOSE       ] = 1;
    iparam[IPARAM_AUTOTUNING    ] = 1;
    iparam[IPARAM_INPUTFMT      ] = 0;
    iparam[IPARAM_OUTPUTFMT     ] = 0;
    iparam[IPARAM_TRACE         ] = 0;
    iparam[IPARAM_DAG           ] = 0;
    iparam[IPARAM_ASYNC         ] = 1;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_RHBLK         ] = 0;
    iparam[IPARAM_TNTPIV_MODE   ] = PLASMA_TOURNAMENT_LU;
    iparam[IPARAM_TNTPIV_SIZE   ] = 4;
    iparam[IPARAM_INPLACE       ] = PLASMA_OUTOFPLACE;

    get_thread_count( &(iparam[IPARAM_THRDNBR]) );

    for (i = 1; i < argc && argv[i]; ++i) {
        if (startswith( argv[i], "--help" )) {
            show_help( argv[0] );
            return EXIT_SUCCESS;
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else if (startswith( argv[i], "--check" )) {
            iparam[IPARAM_CHECK] = 1;
        } else if (startswith( argv[i], "--nocheck" )) {
            iparam[IPARAM_CHECK] = 0;
        } else if (startswith( argv[i], "--warmup" )) {
            iparam[IPARAM_WARMUP] = 1;
        } else if (startswith( argv[i], "--nowarmup" )) {
            iparam[IPARAM_WARMUP] = 0;
        } else if (startswith( argv[i], "--dyn" )) {
            iparam[IPARAM_SCHEDULER] = 1;
        } else if (startswith( argv[i], "--nodyn" )) {
            iparam[IPARAM_SCHEDULER] = 0;
        } else if (startswith( argv[i], "--quark" )) {
            iparam[IPARAM_RUNTIME] = 1;
        } else if (startswith( argv[i], "--ompss" )) {
            iparam[IPARAM_RUNTIME] = 2;
        /* } else if (startswith( argv[i], "--atun" )) { */
        /*     iparam[IPARAM_AUTOTUNING] = 1; */
        /* } else if (startswith( argv[i], "--noatun" )) { */
        /*     iparam[IPARAM_AUTOTUNING] = 0; */
        } else if (startswith( argv[i], "--trace" )) {
            iparam[IPARAM_TRACE] = 1;
        } else if (startswith( argv[i], "--notrace" )) {
            iparam[IPARAM_TRACE] = 0;
        } else if (startswith( argv[i], "--dag" )) {
            iparam[IPARAM_DAG] = 1;
        } else if (startswith( argv[i], "--nodag" )) {
            iparam[IPARAM_DAG] = 0;
        } else if (startswith( argv[i], "--sync" )) {
            iparam[IPARAM_ASYNC] = 0;
        } else if (startswith( argv[i], "--async" )) {
            iparam[IPARAM_ASYNC] = 1;
        } else if (startswith( argv[i], "--n_range=" )) {
            get_range( strchr( argv[i], '=' ) + 1, &start, &stop, &step );
        } else if (startswith( argv[i], "--m=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_M]) );
        } else if (startswith( argv[i], "--n=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_N]) );
            start = iparam[IPARAM_N];
            stop  = iparam[IPARAM_N];
            step  = 1;
        } else if (startswith( argv[i], "--nb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NB]) );
            iparam[IPARAM_MB] = iparam[IPARAM_NB];
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_K]) );
        } else if (startswith( argv[i], "--k=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_K]) );
        } else if (startswith( argv[i], "--ib=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_IB]) );
        } else if (startswith( argv[i], "--ifmt=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_INPUTFMT]) );
        } else if (startswith( argv[i], "--ofmt=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_OUTPUTFMT]) );
        } else if (startswith( argv[i], "--thrdbypb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR_SUBGRP]) );
        } else if (startswith( argv[i], "--niter=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &iparam[IPARAM_NITER] );
        } else if (startswith( argv[i], "--mx=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_MX]) );
        } else if (startswith( argv[i], "--nx=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NX]) );
        } else if (startswith( argv[i], "--rhblk=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_RHBLK]) );
        } else if (startswith( argv[i], "--tntmode=" )) {
            int tmp;
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &tmp );
            iparam[IPARAM_TNTPIV_MODE] = (tmp == 1) ? PLASMA_TOURNAMENT_LU : PLASMA_TOURNAMENT_QR;
        } else if (startswith( argv[i], "--tntsize=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_TNTPIV_SIZE]) );
        } else if (startswith( argv[i], "--inplace" )) {
            iparam[IPARAM_INPLACE] = PLASMA_INPLACE;
        } else if (startswith( argv[i], "--outplace" )) {
            iparam[IPARAM_INPLACE] = PLASMA_OUTOFPLACE;
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }

    m  = iparam[IPARAM_M];
    mx = iparam[IPARAM_MX];
    nx = iparam[IPARAM_NX];

    /* Initialize Plasma */
	if ( iparam[IPARAM_RUNTIME] == PLASMA_QUARK )
		PLASMA_Runtime_Init(omp_get_num_threads(), PLASMA_QUARK);
	else
		PLASMA_Runtime_Init(omp_get_num_threads(), PLASMA_OMPSS);

    if ( iparam[IPARAM_SCHEDULER] ){
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING );
        if ( iparam[IPARAM_RUNTIME] == PLASMA_QUARK )
            PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK );
        else
            PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS );
    }
    else
        PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_STATIC_SCHEDULING );

    /* if ( !iparam[IPARAM_AUTOTUNING] ) { */
    PLASMA_Disable(PLASMA_AUTOTUNING);
    PLASMA_Set(PLASMA_TILE_SIZE,        iparam[IPARAM_NB] );
    PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, iparam[IPARAM_IB] );
    /* } else { */
    /*     PLASMA_Get(PLASMA_TILE_SIZE,        &iparam[IPARAM_NB] ); */
    /*     PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &iparam[IPARAM_IB] ); */
    /* } */

    /* Householder mode */
    if (iparam[IPARAM_RHBLK] < 1) {
        PLASMA_Set(PLASMA_HOUSEHOLDER_MODE, PLASMA_FLAT_HOUSEHOLDER);
    } else {
        PLASMA_Set(PLASMA_HOUSEHOLDER_MODE, PLASMA_TREE_HOUSEHOLDER);
        PLASMA_Set(PLASMA_HOUSEHOLDER_SIZE, iparam[IPARAM_RHBLK]);
    }

    /* Tournament pivoting */
    PLASMA_Set(PLASMA_TNTPIVOTING_MODE, iparam[IPARAM_TNTPIV_MODE]);
    PLASMA_Set(PLASMA_TNTPIVOTING_SIZE, iparam[IPARAM_TNTPIV_SIZE]);

    /* Layout conversion */
    PLASMA_Set(PLASMA_TRANSLATION_MODE, iparam[IPARAM_INPLACE]);

    print_header( argv[0], iparam);

    if (step < 1) step = 1;

    Test( -1, iparam ); /* print header */
    for (i = start; i <= stop; i += step)
    {
        if ( nx > 0 ) {
            iparam[IPARAM_M] = i;
            iparam[IPARAM_N] = max(1, i/nx);
        } else if ( mx > 0 ) {
            iparam[IPARAM_M] = max(1, i/mx);
            iparam[IPARAM_N] = i;
        } else {
            if ( m == -1 )
                iparam[IPARAM_M] = i;
            iparam[IPARAM_N] = i;
        }
        Test( iparam[IPARAM_N], iparam );
    }

    PLASMA_Finalize();

    /* if (gnuplot) { */
    /*         printf( "%s\n%s\n", */
    /*                 "e", */
    /*                 gnuplot > 1 ? "" : "pause 10" ); */
    /* } */

    return EXIT_SUCCESS;
}
