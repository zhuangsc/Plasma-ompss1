      PROGRAM EXAMPLE_SGELS_F
*     
*********************************************************************
*     PLASMA example routine (version 2.6.0)                        
*     Author: Bilel Hadri                                           
*     Release Date: November, 15th 2010                             
*     PLASMA is a software package provided by Univ. of Tennessee,  
*     Univ. of California Berkeley and Univ. of Colorado Denver.    
*     @generated s Tue Jan  7 11:45:21 2014
*********************************************************************
*     
      IMPLICIT NONE
*     
      INCLUDE "plasmaf.h"
*     
*     Purpose
*     =======
*     
*     FORTRAN EXAMPLE FOR PLASMA_SGELS
*     Example for solving a system of linear equations using QR factorization
*     
*     =====================================================================
*     
*     .. Parameters ..
      INTEGER           CORES, M, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( M = 20 )
      PARAMETER         ( N = 15 )
      PARAMETER         ( NRHS = 5 )
      COMPLEX*16        ZONE
      PARAMETER         ( ZONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16  A1( M, N ), B1( MAX(M,N), NRHS )
      COMPLEX*16  A2( M, N ), B2( MAX(M,N), NRHS )
      COMPLEX*16  RISU( MAX(M,N), NRHS)
      DOUBLE PRECISION  RWORK( MAX(M,N ))
      INTEGER           HT( 2 )
      DOUBLE PRECISION  XNORM, ANORM, BNORM, RNORM, EPS
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
*     ..
*     .. External Subroutines ..
      DOUBLE PRECISION  DLAMCH, SLANGE
      EXTERNAL          ZLARNV, DLAMCH, SLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_ALLOC_WORKSPACE_SGELS
      EXTERNAL          PLASMA_SGELS, PLASMA_FINALIZE
      EXTERNAL          PLASMA_DEALLOC_HANDLE
      EXTERNAL          SGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*     
      DO  I = 1, 4
         ISEED( I ) = 1
      ENDDO
*     
*     Initialize Plasma
*     
      CALL PLASMA_INIT( CORES, INFO )
      WRITE(*,*) "-- PLASMA is initialized on", CORES, "cores."
*     
*     Initialization of the matrix
*     
      CALL ZLARNV( 1, ISEED, M*N, A1 )
      A2(:,:)=A1(:,:)
*     
*     Initialization of the RHS
*     
      CALL ZLARNV( 1, ISEED, MAX(M,N)*NRHS, B1 )
      B2(:,:)=B1(:,:)
*     
*     Allocate T
*     
      CALL PLASMA_ALLOC_WORKSPACE_SGELS( M, N, HT, INFO )
*     
*     Perform the QR solve
*     
      CALL PLASMA_SGELS( PlasmaNoTrans, M, N, NRHS,
     &     A2, M, HT, B2, MAX(M,N), INFO )
*     
*     Check the solution
*     
      XNORM = SLANGE('I',MIN(M,N), NRHS, B2, MIN(M,N), RWORK)
      ANORM = SLANGE('I',M, N, A1, M, RWORK)
      BNORM = SLANGE('I',MIN(M,N), NRHS, B1, MIN(M,N), RWORK)

      CALL SGEMM('No transpose','No transpose', M, NRHS, N, ZONE,
     $     A1, M, B2, MAX(M,N), -ZONE, B1, MAX(M,N))

      IF (M >=N ) THEN
         CALL SGEMM('ConjTranspose','No transpose', N, NRHS, M, ZONE,
     $        A1, M, B1, MAX(M,N), -ZONE, RISU, M)
         RNORM = SLANGE('I', M, NRHS, RISU, N, RWORK)
      ELSE
         CALL SGEMM('ConjTranspose','No transpose', N, NRHS, M, ZONE,
     $        A1, M, B1, MAX(M,N), -ZONE, RISU, N)
         RNORM = SLANGE('I', N, NRHS, RISU, N, RWORK)
      ENDIF

      EPS= DLAMCH('Epsilon')

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $     RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      IF ((RNORM > 60.0).AND.( INFO < 0 )) THEN
         WRITE(*,*) "-- Error in SGELS example !"
      ELSE
         WRITE(*,*) "-- Run of SGELS example successful !"
      ENDIF
*     
*     Deallocate T
*     
      CALL PLASMA_DEALLOC_HANDLE( HT, INFO )
*     
*     Finalize Plasma
*     
      CALL PLASMA_FINALIZE( INFO )
*     
*     End of EXAMPLE_SGELS.
*     
      END PROGRAM EXAMPLE_SGELS_F
