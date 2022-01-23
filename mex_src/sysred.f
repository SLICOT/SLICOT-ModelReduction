C SYSRED.F - Gateway function for SLICOT model reduction routines
C            AB09MD.F, AB09ND.F, AB09ED.F, AB09FD.F and AB09GD.F.
C
C RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ar,Br,Cr,Dr,HSV,info] = sysred(meth,A,B,C,D,tol,discr,ord,alpha)
C
C Purpose:
C   To find a reduced order state-space representation (Ar,Br,Cr,Dr)
C   and the Hankel singular values HSV for the alpha-stable part of a
C   continuous- or discrete-time state-space system (A,B,C,D).
C   The order of the reduced model is determined either by the number
C   of Hankel-singular values greater than tol or by the desired
C   order ord.
C
C Input parameters:
C   meth  - method flag with decimal form c*10+m, where:
C           m specifies the basic model reduction method;
C           c specifies the comprime factorization approach to be
C           used in conjunction with the method specified by m.
C           Allowed values for m:
C             m = 1 : Balance & Truncate method with balancing
C             m = 2 : Balance & Truncate method (no balancing)
C             m = 3 : Singular Perturbation Approximation with balancing
C             m = 4 : Singular Perturbation Approximation (no balancing)
C             m = 5 : Optimal Hankel-Norm Approximation.
C           Allowed values for c (only for m = 1..4):
C             c = 0 : no coprime factorization is used (default)
C             c = 1 : RCF with inner denominator
C             c = 2 : LCF with inner denominator
C             c = 3 : RCF with ALPHA stability degree
C             c = 4 : LCF with ALPHA stability degree.
C   A,B,
C   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
C           and P-by-M, respectively.
C   tol   - (optional) tolerance vector for determining the order of
C           reduced system, of the form [tol1, tol2, tol3], where
C             tol1 specifies the tolerance for model reduction;
C                  default: tol1 = epsilon_machine*Hankel_norm(A,B,C)
C             tol2 specifies the tolerance for minimal realization in
C                  case of m = 3, 4 or 5;
C                  default: tol2 = epsilon_machine*Hankel_norm(A,B,C)
C             tol3 specifies the controllability/observability tolerance
C                  for computing coprime factorizations, as follows:
C                  controllability tolerance in case c = 1 or 3;
C                  default: epsilon_machine*max(norm(A),norm(B));
C                  observability tolerance in case c = 2 or 4;
C                  default: epsilon_machine*max(norm(A),norm(C)).
C   discr - (optional) type of system
C              = 0 : continuous-time (default)
C              = 1 : discrete-time.
C   ord   - (optional) desired order of reduced system
C             default: ord = -1 (order determined automatically).
C   alpha - (optional) stability boundary for the eigenvalues of A
C             default:    -sqrt(epsilon_machine)  for continuous-time
C                      1.0-sqrt(epsilon_machine)  for discrete-time.
C
C Output parameters:
C   Ar, Br,
C   Cr, Dr - matrices of the reduced system.
C   HSV    - Hankel singular values of the alpha-stable part.
C   info   - warning message code.
C            info = 1 - selected order greater than the order
C                       of a minimal realization
C            info = 2 - selected order less than the order of
C                       the unstable part.
C
C Contributor:
C   A. Varga, German Aerospace Center,
C   DLR Oberpfaffenhofen, February 1999.
C
C Revisions:
C   V. Sima, February 2001 (bug found by S. Steer); March 2005,
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
      CHARACTER         DICO, EQUIL, FACT, JOBCF, JOBMR, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDWORK, M, N,
     $                  NR, NS, P
      DOUBLE PRECISION  ALPHA, TOL1, TOL2, TOL3
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                DWORK(:), HSV(:)
      INTEGER, ALLOCATABLE::          IWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           COFACT, DISCR
      INTEGER           CFMETH, I, M1, METH, N1, N2, N3, P1
      DOUBLE PRECISION  DUM, TOL(3)
C
C .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          LSAME, DLAMCH
C
C .. External Subroutines ..
      EXTERNAL          AB09MD, AB09ND, AB09ED, AB09FD, AB09GD, DLACPY
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, MOD, SQRT
C
C Check for proper number of arguments.
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSRED requires at least 5 input arguments' )
      ELSE IF( NLHS .GT. 6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SYSRED requires at most 6 output arguments' )
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
         CALL mexErrMsgTxt( 'METH must be an integer scalar 0 or 1' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), DUM, 1 )
      METH = DUM
      CFMETH = METH/10
      METH = MOD( METH, 10 )
      IF( METH.LE.0 .OR. METH.GT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 5 the only admissible values' )
      END IF
      IF( CFMETH.LT.0 .OR. CFMETH.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'Coprime METH has 0 ... 4 the only admissible values' )
      END IF
      COFACT = CFMETH.GT.0
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
C   tol(1x3)
C
      TOL1 = ZERO
      TOL2 = ZERO
      TOL3 = ZERO
      IF( NRHS.GT.5 ) THEN
         I = mxGetM( PRHS(6) )*mxGetN( PRHS(6) )
         IF( I.GT.3 ) THEN
            CALL mexErrMsgTxt
     $           ( 'TOL must be a vector with at most 3 elements' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TOL, I )
         IF( I.GT.0 ) TOL1 = TOL(1)
         IF( I.GT.1 ) TOL2 = TOL(2)
         IF( I.GT.2 ) TOL3 = TOL(3)
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
            IF( ALPHA.EQ.ONE) ALPHA = ONE - SQRT( DLAMCH( 'E' ) )
         ELSE
            IF( ALPHA.EQ.ZERO ) ALPHA = -SQRT( DLAMCH( 'E' ) )
         END IF
      ELSE
         ALPHA = -SQRT( DLAMCH( 'E' ) )
         IF( DISCR ) ALPHA = ONE + ALPHA
      END IF
C
C Determine the lenghts of working arrays.
C
      LWR = MAX( 1, N*( 2*N + MAX( N, M + P ) + 5 ) + ( N*( N+1 ) )/2 )
      IF( COFACT ) THEN
C        Coprime factorization
         IF( CFMETH.EQ.1 ) THEN
C           JOBCF = 'R' and FACT = 'I'
            LDWORK = ( N+M )*( M+P ) + MAX( M*( M+2 ), 4*M, 4*P, LWR )
         ELSE IF( CFMETH.EQ.2 ) THEN
C           JOBCF = 'L' and FACT = 'I'
            LDWORK = N*( 2*MAX( M, P ) + P ) +
     $                    MAX( M, P )*( MAX( M, P ) + P ) +
     $                    MAX( N*P + MAX( N*( N+5 ), P*( P+2 ), 4*P,
     $                    4*M ), LWR )
         ELSE IF( CFMETH.EQ.3 ) THEN
C           JOBCF = 'R' and FACT = 'S'
            LDWORK = ( N+M )*( M+P ) + MAX( 5*M, 4*P, LWR )
         ELSE
C           JOBCF = 'L' and FACT = 'S'
            LDWORK = N*( 2*MAX( M, P ) + P ) +
     $                    MAX( M, P )*( MAX( M, P ) + P ) +
     $                    MAX( N*P + MAX( N*( N+5 ), 5*P, 4*M ), LWR )
         END IF
      ELSE
         IF( METH.EQ.5 ) THEN
C           HNA
            LDWORK = MAX( LWR, N*( M+P+2 ) + 2*M*P + MIN( N, M ) +
     $                        MAX( 3*M + 1, MIN( N, M ) + P ) )
         ELSE
C           B & T or SPA
            LDWORK = LWR
         END IF
      END IF
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
     $           IWORK( MAX( 1, 2*N, M ) ) )
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
      IF( .NOT.COFACT ) THEN
C
C        Balance & Truncate Approximation with balancing.
C
         IF( METH.EQ.1 ) THEN
            CALL AB09MD( DICO, 'Balance', EQUIL, ORDSEL, N, M, P, NR,
     $                   ALPHA, A, LDA, B, LDB, C, LDC, NS, HSV, TOL1,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
C
C        Balance & Truncate Approximation without balancing.
C
         ELSE IF( METH.EQ.2 ) THEN
            CALL AB09MD( DICO, 'No balance', EQUIL, ORDSEL, N, M, P, NR,
     $                   ALPHA, A, LDA, B, LDB, C, LDC, NS, HSV, TOL1,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
C
C        Singular Perturbation Approximation with balancing.
C
         ELSE IF( METH.EQ.3 ) THEN
            CALL AB09ND( DICO, 'Balance', EQUIL, ORDSEL, N, M, P, NR,
     $                   ALPHA, A, LDA, B, LDB, C, LDC, D, LDD, NS, HSV,
     $                   TOL1, TOL2, IWORK, DWORK, LDWORK, IWARN, INFO )
C
C        Singular Perturbation Approximation without balancing.
C
         ELSE IF( METH.EQ.4 ) THEN
            CALL AB09ND( DICO, 'No balance', EQUIL, ORDSEL, N, M, P, NR,
     $                   ALPHA, A, LDA, B, LDB, C, LDC, D, LDD, NS, HSV,
     $                   TOL1, TOL2, IWORK, DWORK, LDWORK, IWARN, INFO )
C
C        Hankel-Norm Approximation.
C
         ELSE IF( METH.EQ.5 ) THEN
            CALL AB09ED( DICO, EQUIL, ORDSEL, N, M, P, NR, ALPHA,
     $                   A, LDA, B, LDB, C, LDC, D, LDD, NS, HSV, TOL1,
     $                   TOL2, IWORK, DWORK, LDWORK, IWARN, INFO )
         END IF
      ELSE
C
C        Coprime factorization approach.
C
         IF( CFMETH.EQ.1 )  THEN
            JOBCF = 'R'
            FACT  = 'I'
         ELSE IF( CFMETH.EQ.2 ) THEN
            JOBCF = 'L'
            FACT  = 'I'
         ELSE IF( CFMETH.EQ.3 ) THEN
            JOBCF = 'R'
            FACT  = 'S'
         ELSE IF( CFMETH.EQ.4 ) THEN
            JOBCF = 'L'
            FACT  = 'S'
         END IF
C
         IF( METH.EQ.1 .OR. METH.EQ.3 ) THEN
            JOBMR = 'B'
         ELSE IF( METH.EQ.2 .OR. METH.EQ.4 ) THEN
            JOBMR = 'N'
         END IF
C
         IF( METH.LE.2 ) THEN
C
C           Combining with Balance & Truncate Approximation.
C
            CALL AB09FD( DICO, JOBCF, FACT, JOBMR, EQUIL, ORDSEL, N, M,
     $                   P, NR, ALPHA, A, LDA, B, LDB, C, LDC, NS, HSV,
     $                   TOL1, TOL3, IWORK, DWORK, LDWORK, IWARN, INFO )
         ELSE IF( METH.LE.4 ) THEN
C
C           Combining with Singular Perturbation Approximation.
C
            CALL AB09GD( DICO, JOBCF, FACT, JOBMR, EQUIL, ORDSEL, N, M,
     $                   P, NR, ALPHA, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   NS, HSV, TOL1, TOL2, TOL3, IWORK, DWORK,
     $                   LDWORK, IWARN, INFO )
         END IF
      END IF
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
     $         CALL dlacpy( 'F', NR, M, B, LDB, B, NR )
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
      DEALLOCATE ( A, B, C, D, HSV, IWORK, DWORK )
C
C Error handling.
C
      IF( INFO.NE.0 ) THEN
         IF( .NOT.COFACT ) THEN
            IF( METH.EQ.1 .OR. METH.EQ.2 ) THEN
               write( TEXT,'( " INFO =", I4, " ON EXIT FROM AB09MD" )' )
     $                INFO
            ELSE IF( METH.EQ.3 .OR. METH.EQ.4 ) THEN
               write( TEXT,'( " INFO =", I4, " ON EXIT FROM AB09ND" )' )
     $                INFO
            ELSE IF( METH.EQ.5 ) THEN
               write( TEXT,'( " INFO =", I4, " ON EXIT FROM AB09ED" )' )
     $                INFO
            END IF
         ELSE
            IF( METH.EQ.1 .OR. METH.EQ.2 ) THEN
               write( TEXT,'( " INFO =", I4, " ON EXIT FROM AB09FD" )' )
     $                INFO
            ELSE IF( METH.EQ.3 .OR. METH.EQ.4 ) THEN
               write( TEXT,'( " INFO =", I4, " ON EXIT FROM AB09GD" )' )
     $                INFO
            END IF
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SYSRED ***
      END
