C CONRED.F - Gateway function for SLICOT controller reduction routine
C            SB16AD.F.
C
C RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Acr,Bcr,Ccr,Dcr,HSVC,info] = CONRED(meth,Ac,Bc,Cc,Dc,A,B,C,D,...
C                                        tol,discr,ord,alpha)
C
C Purpose:
C   To compute a reduced order controller (Acr,Bcr,Ccr,Dcr) for an
C   original state-space controller representation (Ac,Bc,Cc,Dc) by
C   using the frequency-weighted square-root or balancing-free
C   square-root Balance & Truncate (B&T) or Singular Perturbation
C   Approximation (SPA) model reduction methods. The algorithm tries
C   to minimize the norm of the frequency-weighted error
C
C           ||V*(K-Kr)*W||
C
C   where K and Kr are the transfer-function matrices of the original
C   and reduced order controllers, respectively. V and W are special
C   frequency-weighting transfer-function matrices constructed
C   to enforce closed-loop stability and/or closed-loop performance.
C   If G is the transfer-function matrix of the open-loop system, then
C   the following weightings V and W can be used:
C                      -1
C      (a)   V = (I-G*K) *G, W = I - to enforce closed-loop stability;
C                              -1
C      (b)   V = I,  W = (I-G*K) *G - to enforce closed-loop stability;
C                      -1              -1
C      (c)   V = (I-G*K) *G, W = (I-G*K)  - to enforce closed-loop
C            stability and performance.
C
C   G has the state space representation (A,B,C,D).
C   If K is unstable, only the ALPHA-stable part of K is reduced.
C
C Input parameters:
C   meth  - method flag of decimal form ijkl, where
C             i = 1 : use standard choice for controllability Grammian;
C             i = 2 : use stability garanteeing choice for the
C                     controllability Grammian;
C             j = 1 : use standard choice for observability Grammian;
C             j = 2 : use stability garanteeing choice for the
C                     observability Grammian;
C             k = 1 : use the square-root BT method;
C             k = 2 : use the balancing-free square-root BT method;
C             k = 3 : use the square-root SPA method;
C             k = 4 : use the balancing-free square-root SPA method;
C             l = 1 : no weightings are used;
C             l = 2 : stability enforcing left (output) weighting;
C             l = 3 : stability enforcing right (input) weighting;
C             l = 4 : stability and performance enforcing weightings.
C           Note: For a complete explanation on Grammian choices
C                 see subroutine SB16AD.
C   AC,BC,
C   CC,DC - state-space system matrices of size NC-by-NC, NC-by-P,
C           M-by-NC, and M-by-P, respectively.
C   A,B,
C   C,D   - (optional) state-space system matrices of size
C           N-by-N, N-by-M, P-by-N, and P-by-M, respectively.
C   tol   - (optional) tolerance vector for determining the order of
C           reduced system, of the form [tol1, tol2], where:
C             tol1 specifies the tolerance for model reduction.
C                  Default: tol1 = NCS*epsilon_machine*HSVC(1), where
C                  NCS is the order of the alpha-stable part of K.
C             tol2 specifies the tolerance for minimal realization.
C                  Default: tol2 = NCS*epsilon_machine*HSVC(1).
C   discr - (optional) type of system:
C              = 0 : continuous-time (default);
C              = 1 : discrete-time.
C   ord   - (optional) desired order of reduced system.
C              Default: ord = -1 (order determined automatically).
C   alpha - (optional) stability boundary for the eigenvalues of AC.
C              Default:    -sqrt(epsilon_machine)  for continuous-time;
C                       1.0-sqrt(epsilon_machine)  for discrete-time.
C
C Output parameters:
C   Acr, Bcr,
C   Ccr, Dcr - matrices of the reduced controller.
C   HSVC   - frequency-weighted Hankel singular values of the
C            alpha-stable part of K.
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
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
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
      CHARACTER         DICO, EQUIL, JOBC, JOBMR, JOBO, ORDSEL, WEIGHT
      INTEGER           INFO, IWARN, LDA, LDAC, LDB, LDBC, LDC, LDCC,
     $                  LDD, LDDC, LDWORK, M, N, NC, NCR, NCS, P
      DOUBLE PRECISION  ALPHA, TOL1, TOL2
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), AC(:,:), B(:,:), BC(:,:),
     $                                C(:,:), CC(:,:), D(:,:), DC(:,:),
     $                                DWORK(:), HSVC(:)
      INTEGER, ALLOCATABLE::          IWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           DISCR, FRWGHT, SYSTEM
      INTEGER           I, J, K, LW, M1, MC, METH, MP, N1, N2, N3, NNC,
     $                  P1, PC
      DOUBLE PRECISION  DUM, TOL(2)
C
C .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          LSAME, DLAMCH
C
C .. External Subroutines ..
      EXTERNAL          SB16AD, DLACPY
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONRED requires at least 5 input arguments' )
      ELSE IF( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONRED requires at most 6 output arguments' )
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
      I = METH/1000
      METH = METH - I*1000
      J = METH/100
      METH = METH - J*100
      K = METH/10
      METH = METH - K*10
      IF( METH.LT.1 .OR. METH.GT.4 .OR. I.LT.1 .OR. I.GT.2
     $    .OR. J.LT.1 .OR. J.GT.2 .OR. K.LT.1 .OR. K.GT.4 ) THEN
         CALL mexErrMsgTxt
     $   ( 'Invalid value for METH ' )
      END IF
      IF( I.EQ.1 ) THEN
         JOBC = 'S'
      ELSE
         JOBC = 'E'
      END IF
      IF( J.EQ.1 ) THEN
         JOBO = 'S'
      ELSE
         JOBO = 'E'
      END IF
      IF( METH.EQ.1 ) THEN
         WEIGHT = 'N'
      ELSE IF( METH.EQ.2 ) THEN
         WEIGHT = 'O'
      ELSE IF( METH.EQ.3 ) THEN
         WEIGHT = 'I'
      ELSE
         WEIGHT = 'P'
      END IF
      IF( K.EQ.1 ) THEN
         JOBMR = 'B'
      ELSE IF( K.EQ.2 ) THEN
         JOBMR = 'F'
      ELSE IF( K.EQ.3 ) THEN
         JOBMR = 'S'
      ELSE
         JOBMR = 'P'
      END IF
      FRWGHT = METH.GT.1
C
C   AC(NCxNC), BC(NCxP), CC(MxNC), DC(MxP)
C
      NC = mxGetM( PRHS(2) )
      MC = mxGetN( PRHS(3) )
      PC = mxGetM( PRHS(4) )
      N1 = mxGetN( PRHS(2) )
      N2 = mxGetM( PRHS(3) )
      N3 = mxGetN( PRHS(4) )
      P1 = mxGetM( PRHS(5) )
      M1 = mxGetN( PRHS(5) )
C
      IF( N1.NE.NC ) THEN
         CALL mexErrMsgTxt
     $        ( 'AC must be a square matrix' )
      END IF
      IF( N2.NE.NC ) THEN
         CALL mexErrMsgTxt
     $        ( 'BC must have the same row dimension as AC' )
      END IF
      IF( N3.NE.NC ) THEN
         CALL mexErrMsgTxt
     $        ( 'CC must have the same column dimension as AC' )
      END IF
      IF( P1.NE.PC ) THEN
         CALL mexErrMsgTxt
     $        ( 'DC must have the same row dimension as CC' )
      END IF
      IF( M1.NE.MC ) THEN
         CALL mexErrMsgTxt
     $        ( 'DC must have the same column dimension as BC' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'AC must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'BC must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'CC must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'DC must be a real matrix' )
      END IF
      IF( MC.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The controller has no inputs' )
      END IF
      IF( PC.LE.0 ) THEN
         CALL mexErrMsgTxt( 'The controller has no outputs' )
      END IF
C
C   A(NxN), B(NxM), C(PxN), D(PxM)
C
      IF( NRHS.GT.5 .AND. NRHS.LT.9 ) THEN
         CALL mexErrMsgTxt
     $        ( 'CONRED requires at least 9 input arguments' )
      END IF
      N  = mxGetM( PRHS(6) )
      M  = mxGetN( PRHS(7) )
      P  = mxGetM( PRHS(8) )
      N1 = mxGetN( PRHS(6) )
      N2 = mxGetM( PRHS(7) )
      N3 = mxGetN( PRHS(8) )
      P1 = mxGetM( PRHS(9) )
      M1 = mxGetN( PRHS(9) )
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
      SYSTEM = NC*MC*PC.GT.0
      IF( .NOT.SYSTEM .AND. NC+MC+PC.GT.0 ) THEN
         CALL mexErrMsgTxt
     $        ( 'Some dimensions of an empty system are nonzero' )
      END IF
      IF( SYSTEM ) THEN
         IF( MC.NE.P ) THEN
            CALL mexErrMsgTxt
     $        ( 'K must have the same number of inputs as the number'//
     $          ' of outputs of G' )
         END IF
         IF( PC.NE.M ) THEN
            CALL mexErrMsgTxt
     $        ( 'K must have the same number of outputs as the number'//
     $          ' of inputs of G' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'A must be a real matrix' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'B must be a real matrix' )
         END IF
         IF( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'C must be a real matrix' )
         END IF
         IF( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'D must be a real matrix' )
         END IF
      ELSE
         WEIGHT = 'N'
         FRWGHT = .FALSE.
      END IF
C
C   tol(1x2)
C
      TOL1 = ZERO
      TOL2 = ZERO
      IF( NRHS.GT.9 ) THEN
         I = mxGetM( PRHS(10) )*mxGetN( PRHS(10) )
         IF( I.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'TOL must be a vector with at most 2 elements' )
         END IF
         IF( mxIsNumeric( PRHS(10) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(10) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(10) ), TOL, I )
         IF( I.GT.0 ) TOL1 = TOL(1)
         IF( I.GT.1 ) TOL2 = TOL(2)
      END IF
C
C   discr
C
      DICO = 'C'
      IF( NRHS.GT.10 ) THEN
         IF( mxGetM( PRHS(11) ).NE.1 .OR. mxGetN( PRHS(11) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(11) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(11) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'DISCR must be an integer scalar 0 or 1' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(11) ), DUM, 1 )
         IF( DUM.NE.ZERO ) DICO = 'D'
      END IF
      DISCR = LSAME( DICO, 'D' )
C
C   ord
C
      ORDSEL = 'A'
      NCR = NC
      IF( NRHS.GT.11 ) THEN
         IF( mxGetM( PRHS(12) ).NE.1 .OR. mxGetN( PRHS(12) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(12) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(12) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(12) ), DUM, 1 )
         IF( DUM.GE.ZERO ) THEN
            ORDSEL = 'F'
            NCR = DUM
            NCR = MIN( NC, NCR )
         END IF
      END IF
C
C   alpha
C
      IF( NRHS.GT.12 ) THEN
         IF( mxGetM( PRHS(13) ).NE.1 .OR. mxGetN( PRHS(13) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(13) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(13) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(13) ), ALPHA, 1 )
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
C Determine the lenghts of working arrays.
C
      LW  = 1
      NNC = N + NC
      MP  = M + P
      IF( FRWGHT ) THEN
         LW = NNC*( NNC + 2*MP ) +
     $        MAX( NNC*( NNC + MAX( NNC, M, P ) + 7 ), MP*( MP + 4 ) )
      ELSE
         LW = NC*( MAX( M, P ) + 5 )
      END IF
      LW = 2*NC*NC + MAX( 1, LW, NC*( 2*NC + 5 ) )
      LDWORK = LW
C
      LDA = MAX( 1, N )
      LDB = MAX( 1, N )
      LDC = MAX( 1, P )
      LDD = MAX( 1, P )
      LDAC = MAX( 1, NC )
      LDBC = MAX( 1, NC )
      LDCC = MAX( 1, PC )
      LDDC = MAX( 1, PC )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
      ALLOCATE ( A( LDA, MAX( 1, N ) ), B( LDB, MAX( 1, M ) ),
     $           C( LDC, MAX( 1, N ) ), D( LDD, MAX( 1, M ) ),
     $           AC( LDAC, MAX( 1, NC ) ), BC( LDBC, MAX( 1, MC ) ),
     $           CC( LDCC, MAX( 1, NC ) ), DC( LDDC, MAX( 1, MC ) ),
     $           DWORK( LDWORK ), HSVC( MAX( 1, NC ) ),
     $           IWORK( MAX( 1, 2*NC, 2*MP ) ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), AC, NC*NC )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), BC, NC*P )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), CC, M*NC )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), DC, M*P )
      IF( NRHS.GT.5 .AND. SYSTEM ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), A, N*N )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), B, N*M )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), C, P*N )
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), D, P*M )
      END IF
C
C Do the actual computations.
C
      EQUIL = 'S'
      CALL SB16AD( DICO, JOBC, JOBO, JOBMR, WEIGHT, EQUIL, ORDSEL,
     $             N, M, P, NC, NCR, ALPHA, A, LDA, B, LDB,
     $             C, LDC, D, LDD, AC, LDAC, BC, LDBC, CC, LDCC,
     $             DC, LDDC, NCS, HSVC, TOL1, TOL2, IWORK, DWORK,
     $             LDWORK, IWARN, INFO )
C
C Copy output to MATLAB workspace.
C
      IF( INFO.EQ.0 ) THEN
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( NCR, NCR, 0 )
            IF( NCR.LT.NC .AND. NCR.GT.0 )
     $         CALL DLACPY( 'F', NCR, NCR, AC, LDAC, AC, NCR )
            CALL mxCopyReal8ToPtr( AC, mxGetPr( PLHS(1) ), NCR*NCR )
         END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NCR, P, 0 )
            IF( NCR.LT.NC .AND. NCR.GT.0 )
     $         CALL DLACPY( 'F', NCR, P, BC, LDBC, BC, NCR )
            CALL mxCopyReal8ToPtr( BC, mxGetPr( PLHS(2) ), NCR*P )
         END IF
         IF( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( M, NCR, 0 )
            CALL mxCopyReal8ToPtr( CC, mxGetPr( PLHS(3) ), M*NCR )
         END IF
         IF( NLHS.GE.4 ) THEN
            PLHS(4) = mxCreateDoubleMatrix( M, P, 0 )
            CALL mxCopyReal8ToPtr( DC, mxGetPr( PLHS(4) ), M*P )
         END IF
C
         IF( NLHS.GE.5 ) THEN
            PLHS(5) = mxCreateDoubleMatrix( NCS, 1, 0 )
            CALL mxCopyReal8ToPtr( HSVC, mxGetPr( PLHS(5) ), NCS )
         END IF
         IF( NLHS.GE.6 ) THEN
            PLHS(6) = mxCreateDoubleMatrix( 1, 1, 0 )
            HSVC(1) = IWARN
            CALL mxCopyReal8ToPtr( HSVC, mxGetPr( PLHS(6) ), 1 )
         END IF
C
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, D, AC, BC, CC, DC, HSVC, IWORK, DWORK )
C
C Error handling.
C
      IF( INFO.NE.0 ) THEN
         WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM SB16AD" )' ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of CONRED ***
      END
