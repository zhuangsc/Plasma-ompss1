/*
    -- MAGMA (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlacpy.cu normal z -> d, Fri Jan 30 19:00:09 2015
       @author Mark Gates
       @author Azzam Haidar
*/

#include "common.h"
#include "runtime.h"                                                                     
#include "core_blas-gpu.h"                                                                     

#define BLK_X 64
#define BLK_Y 32

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to dlaset.
*/
static __device__
void dlacpy_full_device(
    int m, int n,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    /* do only rows inside matrix */
    if ( ind < m ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
    }
}


/*
    Similar to dlacpy_full, but updates only the diagonal and below.
    Blocks that are fully above the diagonal exit immediately.

    Code similar to dlaset.
*/
static __device__
void dlacpy_lower_device(
    int m, int n,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (below diag) */
    bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y));
    /* do only rows inside matrix, and blocks not above diag */
    if ( ind < m && ind + BLK_X > iby ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n && ind >= iby+j; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
    }
}


/*
    Similar to dlacpy_full, but updates only the diagonal and above.
    Blocks that are fully below the diagonal exit immediately.

    Code similar to dlaset.
*/
static __device__
void dlacpy_upper_device(
    int m, int n,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (above diag) */
    bool full = (iby + BLK_Y <= n && (ind + BLK_X <= iby));
    /* do only rows inside matrix, and blocks not below diag */
    if ( ind < m && ind < iby + BLK_Y ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( ind <= iby+j ) {
                    dB[j*lddb] = dA[j*ldda];
                }
            }
        }
    }
}

/*
    kernel wrapper to call the device function.
*/
__global__
void dlacpy_full_kernel(
    int m, int n,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    dlacpy_full_device(m, n, dA, ldda, dB, lddb);
}

__global__
void dlacpy_lower_kernel(
    int m, int n,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    dlacpy_lower_device(m, n, dA, ldda, dB, lddb);
}

__global__
void dlacpy_upper_kernel(
    int m, int n,
    const double *dA, int ldda,
    double       *dB, int lddb )
{
    dlacpy_upper_device(m, n, dA, ldda, dB, lddb);
}


/**
    Purpose
    -------
    DLACPY_Q copies all or part of a two-dimensional matrix dA to another
    matrix dB.
    
    This is the same as DLACPY, but adds queue argument.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_uplo_t
            Specifies the part of the matrix dA to be copied to dB.
      -     = MagmaUpper:      Upper triangular part
      -     = MagmaLower:      Lower triangular part
            Otherwise:  All of the matrix dA
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    dA      DOUBLE_PRECISION array, dimension (LDDA,N)
            The m by n matrix dA.
            If UPLO = MagmaUpper, only the upper triangle or trapezoid is accessed;
            if UPLO = MagmaLower, only the lower triangle or trapezoid is accessed.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[out]
    dB      DOUBLE_PRECISION array, dimension (LDDB,N)
            The m by n matrix dB.
            On exit, dB = dA in the locations specified by UPLO.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_daux2
    ********************************************************************/
int cublasDlacpy_q( PLASMA_enum uplo, int m, int n, double *dA, int ldda, double *dB, int lddb, cudaStream_t stream )
{
    int info = 0;
    /*
    */
    if ((uplo != PlasmaUpperLower) && (uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        info = -1;
        coreblas_error(  -(info), "uplo is invalid" );
        return info;
    }
    if ( m < 0 ){
        info = -2;
        coreblas_error(  -(info), "M should be positive" );
        return info;
    }
    else if ( n < 0 ){
        info = -3;
        coreblas_error(  -(info), "N should be positive" );
        return info;
    }
    else if ( ldda < max(1,m)){
        info = -5;
        coreblas_error(  -(info), "LDA is invalid" );
        return info;
    }
    else if ( lddb < max(1,m)){
        info = -7;
        coreblas_error(  -(info), "LDB is invalid" );
        return info;
    }
    
    if ( info != 0 ) {
        return -1;
    }
    
    if ( m == 0 || n == 0 )
        return 1;
    
    dim3 threads( BLK_X, 1 );
    dim3 grid( (m + BLK_X - 1)/BLK_X, (n + BLK_Y - 1)/BLK_Y );
    /* 
    */ 
    if ( uplo == PlasmaLower ) {
        dlacpy_lower_kernel<<< grid, threads, 0, stream >>> ( m, n, dA, ldda, dB, lddb );
    }
    else if ( uplo == PlasmaUpper ) {
        dlacpy_upper_kernel<<< grid, threads, 0, stream >>> ( m, n, dA, ldda, dB, lddb );
    }
    else {
        dlacpy_full_kernel <<< grid, threads, 0, stream >>> ( m, n, dA, ldda, dB, lddb );
    }
}

/**
    @see magmablas_dlacpy_q
    @ingroup magma_daux2
    ********************************************************************/
int cublasDlacpy(
    cudaStream_t stream, PLASMA_enum uplo, int m, int n,
    double *dA, int ldda,
    double *dB, int lddb )
{
    cublasDlacpy_q( uplo, m, n, dA, ldda, dB, lddb, stream );
}
