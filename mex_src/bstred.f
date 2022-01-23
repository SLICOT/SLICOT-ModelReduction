C BSTRED.F - Gateway function for SLICOT model reduction routine
C            AB09HD.F
C
C RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ar,Br,Cr,Dr,HSV,info] = BSTRED(meth,A,B,C,D,tol,discr,ord,alpha,beta)
C
C Purpose:
C   To find a reduced order state-space system Gr = (Ar,Br,Cr,Dr)
C   from a continuous- or discrete-time original system G = (A,B,C,D)
C   using the balanced stochastic truncation (BST) or the balanced
C   stochastic singular perturbation approximation (BS-SPA) methods.
C   The order of the reduced model is determined either by the number
C   of stochastic Hankel-singular values HSV greater than tol or
C   by the desired order ord.
C
C Input parameters:
C   meth  - method flag to specify the basic model reduction method;
C           Allowed values for meth are:
C             meth = 1 : BST method with balancing;
C             meth = 2 : BST method (no balancing);
C             meth = 3 : BS-SPA method with balancing;
C             meth = 4 : BS-SPA (no balancing).
C   A,B,
C   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
C           and P-by-M, respectively.
C   tol   - (optional) tolerance vector for determining the order of
C           reduced system, of the form [tol1, tol2], where:
C             tol1 specifies the tolerance for model reduction.
C                  Default: tol1 = NS*epsilon_machine, where NS is the
C                  order of the alpha-stable part of G.
C             tol2 specifies the tolerance for computing a minimal
C                  realization when meth = 3 or 4.
C                  Default: tol2 = NS*epsilon_machine.
C   discr - (optional) type of system:
C              = 0 : continuous-time (default);
C              = 1 : discrete-time.
C   ord   - (optional) desired order of reduced system.
C             Default: ord = -1 (order determined automatically).
C   alpha - (optional) stability boundary for the eigenvalues of A.
C             Default:    -sqrt(epsilon_machine)  for continuous-time;
C                      1.0-sqrt(epsilon_machine)  for discrete-time.
C   beta  - (optional) absolute/relative error weighting parameter.
C           beta must be positive if D has not a full row rank.
C             Default: 0 (pure relative method).
C
C Output parameters:
C   Ar, Br,
C   Cr, Dr - matrices of the reduced system.
C   HSV    - Hankel singular values of the alpha-stable part.
C   info   - warning message code:
C            info = 1 - selected order greater than the order
C                       of a minimal realization;
C            info = 2 - selected order corresponds to repeated singular
C                       values, which are neither all included nor all
C                       excluded from the reduced model;
C            info = 3 - selected order less than the order of
C                       the unstable part.
C
C Contributors:
C   D. Sima, University of Bucharest, and
C   A. Varga, German Aerospace Center,
C   DLR Oberpfaffenhofen, March 2001.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C   Apr. 2009, Dec. 2012.
C
C  *********************************************************************
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, TWO, ZERO
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0 )
C
C .. Mex-file interface parameters ..
      INTEGER           PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      INTEGER           mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         DICO, EQUIL, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDWORK, M, N,
     $                  NR, NS, P
      DOUBLE PRECISION  ALPHA, BETA, TOL1, TOL2
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                DWORK(:), HSV(:)
      INTEGER, ALLOCATABLE::          IWORK(:)
      LOGICAL, ALLOCATABLE::          BWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           DISCR
      INTEGER           I, M1, MB, METH, N1, N2, N3, P1
      DOUBLE PRECISION  DUM, TOL(2)
C
C .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          LSAME, DLAMCH
C
C .. External Subroutines ..
      EXTERNAL          AB09HD, DLACPY
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'BSTRED requires at least 5 input arguments' )
      ELSE IF( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'BSTRED requires at most 6 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   meth
C
      IF( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be a scalar' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), DUM, 1 )
      METH = DUM
      IF( METH.LE.0 .OR. METH.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 4 the only admissible values' )
      END IF
      IF( METH.EQ.1 ) THEN
         JOB = 'B'
      ELSE IF( METH.EQ.2 ) THEN
         JOB = 'F'
      ELSE IF( METH.EQ.3 ) THEN
         JOB = 'S'
      ELSE
         JOB = 'P'
      END IF
C
C   A(NxN), B(NxM), C(PxN), D(PxM)
C
      N  = mxGetM( PRHS(2) )
      M  = mxGetN( PRHS(3) )
      P  = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(2) )
      N2 = mxGetM( PRHS(3) )
      N3 = mxGetN( PRHS(4) )
      P1 = mxGetM( PRHS(5) )
      M1 = mxGetN( PRHS(5) )
C
      IF( N1.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'A must be a square matrix' )
      END IF
      IF( N2.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A' )
      END IF
      IF( N3.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A' )
      END IF
      IF( P1.NE.P ) THEN
         CALL mexErrMsgTxt
     $        ( 'D must have the same row dimension as C' )
      END IF
      IF( M1.NE.M ) THEN
         CALL mexErrMsgTxt
     $        ( 'D must have the same column dimension as B' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
      IF( M.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no inputs' )
      END IF
      IF( P.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The system has no outputs' )
      END IF
C
C   tol(1x2)
C
      TOL1 = ZERO
      TOL2 = ZERO
      IF( NRHS.GT.5 ) THEN
         I = mxGetM( PRHS(6) )*mxGetN( PRHS(6) )
         IF( I.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'TOL must be a vector with at most 2 elements' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TOL, I )
         IF( I.GT.0 ) TOL1 = TOL(1)
         IF( I.GT.1 ) TOL2 = TOL(2)
      END IF
C
C   discr
C
      DICO = 'C'
      IF( NRHS.GT.6 ) THEN
         IF( mxGetM( PRHS(7) ).NE.1 .OR. mxGetN( PRHS(7) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'DISCR must be an integer scalar 0 or 1' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), DUM, 1 )
         IF( DUM.NE.ZERO ) DICO = 'D'
      END IF
      DISCR = LSAME( DICO, 'D' )
C
C   ord
C
      ORDSEL = 'A'
      NR = 0
      IF( NRHS.GT.7 ) THEN
         IF( mxGetM( PRHS(8) ).NE.1 .OR. mxGetN( PRHS(8) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), DUM, 1 )
         IF( DUM.GE.ZERO ) THEN
            ORDSEL = 'F'
            NR = DUM
            NR = MIN( N, NR )
         END IF
      END IF
C
C   alpha
C
      IF( NRHS.GT.8 ) THEN
         IF( mxGetM( PRHS(9) ).NE.1 .OR. mxGetN( PRHS(9) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), ALPHA, 1 )
         IF( DISCR ) THEN
            IF( ALPHA.EQ.ONE )  ALPHA = ONE - SQRT( DLAMCH( 'E' ) )
         ELSE
            IF( ALPHA.EQ.ZERO ) ALPHA = -SQRT( DLAMCH( 'E' ) )
         END IF
      ELSE
         ALPHA = -SQRT( DLAMCH( 'E' ) )
         IF( DISCR ) ALPHA = ONE + ALPHA
      END IF
C
C   beta
C
      BETA = ZERO
      IF( NRHS.GT.9 ) THEN
         IF( mxGetM( PRHS(10) ).NE.1 .OR. mxGetN( PRHS(10) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'BETA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(10) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(10) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'BETA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(10) ), BETA, 1 )
      END IF
C
C Determine the lenghts of working arrays.
C
      MB = M
      IF( BETA.NE.ZERO ) MB = M + P
      LDWORK = 2*N*N + MB*(N+P) + MAX( 2, N*(MAX( N, MB, P )+5),
     $                 2*N*P + MAX( P*(MB+2), 10*N*(N+1) ) )
C
      LDA = MAX( 1, N )
      LDB = MAX( 1, N )
      LDC = MAX( 1, P )
      LDD = MAX( 1, P )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
      ALLOCATE ( A( LDA, MAX( 1, N ) ), B( LDB, MAX( 1, M ) ),
     $           C( LDC, MAX( 1, N ) ), D( LDD, MAX( 1, M ) ),
     $           DWORK( LDWORK ), HSV( MAX( 1, N ) ),
     $           IWORK( MAX( 1, 2*N ) ), BWORK( MAX( 1, 2*N ) ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), D, P*M )
C
C Do the actual computations.
C
      EQUIL = 'S'
      CALL AB09HD( DICO, JOB, EQUIL, ORDSEL, N, M, P, NR,
     $             ALPHA, BETA, A, LDA, B, LDB, C, LDC, D, LDD,
     $             NS, HSV, TOL1, TOL2, IWORK, DWORK, LDWORK,
     $             BWORK, IWARN, INFO )
C
C Copy output to MATLAB workspace.
C
      IF( INFO.EQ.0 ) THEN
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( NR, NR, 0 )
            IF( NR.LT.N .AND. NR.GT.0 )
     $         CALL DLACPY( 'F', NR, NR, A, LDA, A, NR )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), NR*NR )
         END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NR, M, 0 )
            IF( NR.LT.N .AND. NR.GT.0 )
     $         CALL DLACPY( 'F', NR, M, B, LDB, B, NR )
            CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(2) ), NR*M )
         END IF
         IF( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( P, NR, 0 )
            CALL mxCopyReal8ToPtr( C, mxGetPr( PLHS(3) ), P*NR )
         END IF
         IF( NLHS.GE.4 ) THEN
            PLHS(4) = mxCreateDoubleMatrix( P, M, 0 )
            CALL mxCopyReal8ToPtr( D, mxGetPr( PLHS(4) ), P*M )
         END IF
C
         IF( NLHS.GE.5 ) THEN
            PLHS(5) = mxCreateDoubleMatrix( NS, 1, 0 )
            CALL mxCopyReal8ToPtr( HSV, mxGetPr( PLHS(5) ), NS )
         END IF
         IF( NLHS.GE.6 ) THEN
            PLHS(6) = mxCreateDoubleMatrix( 1, 1, 0 )
            HSV(1) = IWARN
            CALL mxCopyReal8ToPtr( HSV, mxGetPr( PLHS(6) ), 1 )
         END IF
C
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, D, HSV, IWORK, DWORK, BWORK )
C
C Error handling.
C
      IF( INFO.NE.0 ) THEN
         WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM AB09HD" )' ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of BSTRED ***
      END
