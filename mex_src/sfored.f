C SFORED.F - Gateway function for SLICOT controller reduction routines
C            SB16BD.F and SB16CD.F.
C
C RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ac,Bc,Cc,Dc,HSV,info] = SFORED(meth,A,B,C,D,F,G,tol,discr,ord)
C
C Purpose:
C   To compute, for a given open-loop model (A,B,C,D), and for
C   given state feedback gain F and full observer gain G,
C   such that A+B*F and A+G*C are stable, a reduced order
C   controller model (Ac,Bc,Cc,Dc) using a coprime factorization
C   based controller reduction approach. For reduction of
C   coprime factors, optionally a stability enforcing frequency-weighted
C   model reduction can be used.  For reduction, either the square-root
C   or the balancing-free square-root versions of the Balance & Truncate
C   (B&T) or Singular Perturbation Approximation (SPA)
C   (only in the non-weighted case) model reduction methods are used in
C   conjunction with stable coprime factorization techniques.
C   The order of the reduced model is determined either by the number
C   of Hankel-singular values greater than tol or by the desired
C   order ord.
C
C Input parameters:
C   meth  - method flag of decimal form mfr to specify
C           the reduction method. The allowed values for m, f and r are:
C             m = 1 : standard coprime factorization;
C             m = 2 : coprime factorization with frequency-weighting;
C             f = 1 : use left coprime factorization;
C             f = 2 : use right coprime factorization;
C             r = 1 : B&T method with balancing;
C             r = 2 : B&T method (no balancing);
C             r = 3 : SPA method with balancing (only for m = 1);
C             r = 4 : SPA (no balancing) (only for m = 1).
C   A,B,
C   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
C           and P-by-M, respectively.
C   F,G   - state-feedack and output-injection matrices of size M-by-N
C           and N-by-P, respectively.
C   tol   - (optional) tolerance vector for determining the order of
C           reduced system, of the form [tol1, tol2], where:
C             tol1 specifies the tolerance for model reduction.
C                  Default: tol1 = N*epsilon_machine*HSV(1),
C                  where HSV(1) is the largest Hankel-singular value
C                  of the extended system Ge (see SB16BD) for m = 1 or
C                  the largest frequecy-weighted Hankel-singular value
C                  (see SB16CD) for m = 2.
C             tol2 specifies, for m = 1, the tolerance for minimal
C                  realization.
C                  Default: tol2 = N*epsilon_machine*HSV(1).
C   discr - (optional) type of system:
C             = 0 : continuous-time (default);
C             = 1 : discrete-time.
C   ord   - (optional) desired order of reduced system.
C             Default: ord = -1 (order determined automatically).
C
C Output parameters:
C   Ac, Bc,
C   Cc, Dc - matrices of the reduced order controller.
C   HSV    - Hankel-singular values of the extended system Ge, for
C            m = 1 (see SB16BD), or the frequency-weighted
C            Hankel-singular values, for m = 2 (see SB16CD).
C   info   - warning message code:
C            info = 1 - selected order greater than the order
C                       of a minimal realization of the controller;
C            info = 2 - selected order corresponds to repeated singular
C                       values, which are neither all included nor all
C                       excluded from the reduced model.
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
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
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
      CHARACTER         DICO, EQUIL, JOBCF, JOBD, JOBMR, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDDC, LDF, LDG,
     $                  LDWORK, M, N, NCR, P
      DOUBLE PRECISION  TOL1, TOL2
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:),  B(:,:),   C(:,:), D(:,:),
     $                                DC(:,:), DWORK(:), F(:,:), G(:,:),
     $                                HSV(:)
      INTEGER, ALLOCATABLE::          IWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           DISCR, FRWGHT, LEFT
      INTEGER           I, IPN, J, LIWORK, LW, M1, M2, METH, N1, N2, N3,
     $                  NF, NG, P1, P2
      DOUBLE PRECISION  DUM, TOL(2)
C
C .. External functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C
C .. External Subroutines ..
      EXTERNAL          SB16BD, SB16CD, DLACPY, DLASET
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
C Check for proper number of arguments.
C
      IF( NRHS.LT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SFORED requires at least 7 input arguments' )
      ELSE IF( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SFORED requires at most 6 output arguments' )
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
      I = METH/100
      METH = METH - I*100
      J = METH/10
      METH = METH - J*10
      IF( METH.LT.1 .OR. ( I.EQ.1 .AND. METH.GT.4 )
     $              .OR. ( I.EQ.2 .AND. METH.GT.2 ) .OR.
     $       I.LT.1 .OR. I.GT.2 .OR. J.LT.1 .OR. J.GT.2 ) THEN
         CALL mexErrMsgTxt('Invalid value for METH ' )
      END IF
      IF( I.EQ.1 ) THEN
         FRWGHT = .FALSE.
      ELSE
         FRWGHT = .TRUE.
      END IF
      IF( J.EQ.1 ) THEN
         JOBCF = 'L'
         LEFT = .TRUE.
      ELSE
         JOBCF = 'R'
         LEFT = .FALSE.
      END IF
      IF( METH.EQ.1 ) THEN
         JOBMR = 'B'
      ELSE IF( METH.EQ.2 ) THEN
         JOBMR = 'F'
      ELSE IF( METH.EQ.3 ) THEN
         JOBMR = 'S'
      ELSE
         JOBMR = 'P'
      END IF
C
C   A(NxN), B(NxM), C(PxN), D(PxM), F(MxN), G(NxP)
C
      N  = mxGetM( PRHS(2) )
      M  = mxGetN( PRHS(3) )
      P  = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(2) )
      N2 = mxGetM( PRHS(3) )
      N3 = mxGetN( PRHS(4) )
      P1 = mxGetM( PRHS(5) )
      M1 = mxGetN( PRHS(5) )
      M2 = mxGetM( PRHS(6) )
      NF = mxGetN( PRHS(6) )
      NG = mxGetM( PRHS(7) )
      P2 = mxGetN( PRHS(7) )
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
      IF( M2.NE.M ) THEN
         CALL mexErrMsgTxt
     $        ( 'F must have the same row dimension as the column '//
     $          'dimension of B' )
      END IF
      IF( NF.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'F must have the same column dimension as A' )
      END IF
      IF( P2.NE.P ) THEN
         CALL mexErrMsgTxt
     $        ( 'G must have the same column dimension as the row '//
     $          'dimension of C' )
      END IF
      IF( NG.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'G must have the same row dimension as A' )
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
      IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'F must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(7) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'G must be a real matrix' )
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
      IPN  = 7
      TOL1 = ZERO
      TOL2 = ZERO
      IF( NRHS.GT.IPN ) THEN
         IPN = IPN + 1
         I = mxGetM( PRHS(IPN) )*mxGetN( PRHS(IPN) )
         IF( I.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'TOL must be a vector with at most 2 elements' )
         END IF
         IF( mxIsNumeric( PRHS(IPN) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(IPN) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPN) ), TOL, I )
         IF( I.GT.0 ) TOL1 = TOL(1)
         IF( I.GT.1 ) TOL2 = TOL(2)
      END IF
C
C   discr
C
      DICO = 'C'
      IF( NRHS.GT.IPN ) THEN
         IPN = IPN + 1
         I = mxGetM( PRHS(IPN) )*mxGetN( PRHS(IPN) )
         IF( I.GT.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(IPN) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(IPN) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'DISCR must be an integer scalar 0 or 1' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPN) ), DUM, 1 )
         IF( DUM.NE.ZERO ) DICO = 'D'
      END IF
      DISCR = LSAME( DICO, 'D' )
C
C   ord
C
      ORDSEL = 'A'
      NCR = 0
      IF( NRHS.GT.IPN ) THEN
         IPN = IPN + 1
         I = mxGetM( PRHS(IPN) )*mxGetN( PRHS(IPN) )
         IF( I.GT.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(IPN) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(IPN) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IPN) ), DUM, 1 )
         IF( DUM.GE.ZERO ) THEN
            ORDSEL = 'F'
            NCR = DUM
            NCR = MIN( N, NCR )
         END IF
      END IF
C
C Determine the lenghts of working arrays.
C
      IF( .NOT.FRWGHT ) THEN
         LW  = MAX( 1, N*( 2*N + MAX( N, M+P ) + 5 ) + ( N*(N+1) )/2 )
         IF( LEFT ) THEN
            LDWORK = ( N + M )*( M + P ) + MAX( LW, 4*M )
            LIWORK = MAX( 1, 2*N, M )
         ELSE
            LDWORK = ( N + P )*( M + P ) + MAX( LW, 4*P )
            LIWORK = MAX( 1, 2*N, P )
         END IF
      ELSE
         IF( LEFT ) THEN
            I = M
         ELSE
            I = P
         END IF
         LDWORK = 2*N*N + MAX( 1, 2*N*N + 5*N, N*MAX( M, P ),
     $                   N*( N + MAX( N, I ) + MIN( N, I ) + 6 ) )
         LIWORK = MAX( 1, N )
      END IF
C
      LDA  = MAX( 1, N )
      LDB  = LDA
      LDC  = MAX( 1, P )
      LDD  = LDC
      LDF  = MAX( 1, M )
      LDG  = LDA
      LDDC = LDF
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
      ALLOCATE ( A( LDA, MAX( 1, N ) ), B( LDB, MAX( 1, M ) ),
     $           C( LDC, MAX( 1, N ) ), D( LDD, MAX( 1, M ) ),
     $           F( LDF, MAX( 1, N ) ), G( LDG, MAX( 1, P ) ),
     $           DC( LDDC, MAX( 1, P ) ), DWORK( LDWORK ),
     $           HSV( MAX( 1, N ) ), IWORK( LIWORK ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), D, P*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), F, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), G, P*N )
C
C Do the actual computations.
C
      JOBD = 'D'
      IF( .NOT.FRWGHT) THEN
         EQUIL = 'S'
         CALL SB16BD( DICO, JOBD, JOBMR, JOBCF, EQUIL, ORDSEL, N,
     $                M, P, NCR, A, LDA, B, LDB, C, LDC, D, LDD,
     $                F, LDF, G, LDG, DC, LDDC, HSV, TOL1, TOL2,
     $                IWORK, DWORK, LDWORK, IWARN, INFO )
      ELSE
         CALL SB16CD( DICO, JOBD, JOBMR, JOBCF, ORDSEL, N, M, P,
     $                NCR, A, LDA, B, LDB, C, LDC, D, LDD, F, LDF,
     $                G, LDG, HSV, TOL1, IWORK, DWORK, LDWORK,
     $                IWARN, INFO )
         CALL DLASET( 'Full', M, P, ZERO, ZERO, DC, LDDC )
      END IF
C
C Copy output to MATLAB workspace.
C
      IF( INFO.EQ.0 ) THEN
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( NCR, NCR, 0 )
            IF( NCR.LT.N .AND. NCR.GT.0 )
     $         CALL DLACPY( 'F', NCR, NCR, A, LDA, A, NCR )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), NCR*NCR )
         END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NCR, P, 0 )
            IF( NCR.LT.N .AND. NCR.GT.0 )
     $         CALL DLACPY( 'F', NCR, P, G, LDG, G, NCR )
            CALL mxCopyReal8ToPtr( G, mxGetPr( PLHS(2) ), NCR*P )
         END IF
         IF( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( M, NCR, 0 )
            CALL mxCopyReal8ToPtr( F, mxGetPr( PLHS(3) ), M*NCR )
         END IF
         IF( NLHS.GE.4 ) THEN
            PLHS(4) = mxCreateDoubleMatrix( M, P, 0 )
            CALL mxCopyReal8ToPtr( DC, mxGetPr( PLHS(4) ), M*P )
         END IF
C
         IF( NLHS.GE.5 ) THEN
            PLHS(5) = mxCreateDoubleMatrix( N, 1, 0 )
            CALL mxCopyReal8ToPtr( HSV, mxGetPr( PLHS(5) ), N )
         END IF
         IF( NLHS.GE.6 ) THEN
            PLHS(6) = mxCreateDoubleMatrix( 1, 1, 0 )
            HSV(1)  = IWARN
            CALL mxCopyReal8ToPtr( HSV, mxGetPr( PLHS(6) ), 1 )
         END IF
C
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, D, F, G, DC, HSV, IWORK, DWORK )
C
C Error handling.
C
      IF( INFO.NE.0 ) THEN
         IF( FRWGHT ) THEN
            WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM SB16CD" )' )
     $                INFO
         ELSE
            WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM SB16BD" )' )
     $                INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SFORED ***
      END
