C FWERED.F - Gateway function for SLICOT model reduction routine
C            AB09ID.F.
C
C RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ar,Br,Cr,Dr,HSV,info] = FWERED(meth,A,B,C,D,V,W,tol,discr,ord,alpha)
C
C Purpose:
C   To find for a given state-space system (A,B,C,D) with
C   transfer-function matrix G a reduced order state-space
C   representation (Ar,Br,Cr,Dr) with the transfer-function matrix Gr
C   such that the frequency-weighted approximation error V*(G-Gr)*W is
C   minimized, where V and W are stable frequency-weighting
C   transfer-function matrices.
C   The frequency-weighted Hankel singular values HSV are also computed.
C   If G is unstable, the approximation is only performed for the
C   alpha-stable part of G.
C   The order of the reduced model is determined either by the number
C   of Hankel-singular values greater than tol or by the desired
C   order ord.
C   A frequency-weighted Balance & Truncate (BT) or
C   Singular Perturbation Approximation (SPA) approach is used.
C
C Input parameters:
C   meth  - method flag of decimal form ijk, where
C             i = 1 : use standard choice for controllability Grammian;
C             i = 2 : use stability garanteeing choice for the
C                     controllability Grammian;
C             j = 1 : use standard choice for observability Grammian;
C             j = 2 : use stability garanteeing choice for the
C                     observability Grammian;
C             k = 1 : use the square-root BT method;
C             k = 2 : use the balancing-free square-root BT method;
C             k = 3 : use the square-root SPA method;
C             k = 4 : use the balancing-free square-root BT method.
C           Note: For a complete explanation on Grammian choices
C                 see subroutine AB09ID.
C   A,B,
C   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
C           and P-by-M, respectively.
C   V     - (optional) contains the system matrix [AV BV; CV DV] of
C           the left weighting; default: V = I.
C   W     - (optional) contains the system matrix [AW BW; CW DW] of
C           the right weighting; default: W = I.
C   tol   - (optional) tolerance vector for determining the order of
C           reduced system, of the form [tol1, tol2], where:
C             tol1 specifies the tolerance for model reduction.
C                  Default: tol1 = NS*epsilon_machine*HSV(1), where NS
C                  is the order of the alpha-stable part of G.
C             tol2 specifies the tolerance for minimal realization.
C                  Default: tol2 = NS*epsilon_machine*HSV(1).
C   discr - (optional) type of system:
C             = 0 : continuous-time (default);
C             = 1 : discrete-time.
C   ord   - (optional) desired order of reduced system.
C             Default: ord = -1 (order determined automatically).
C   alpha - (optional) vector of the form [sdeg, alphac, alphao], where
C             sdeg is the stability boundary for the eigenvalues of A;
C               default:    -sqrt(epsilon_machine)  for continuous-time;
C                        1.0-sqrt(epsilon_machine)  for discrete-time;
C             alphac is the weighting factor for the frequency-weighted
C                    controllability Grammian; default: alphac = 0;
C             alphao is the weighting factor for the frequency-weighted
C                    observability Grammian;   default: alphao = 0.
C             Note: For further explanation of alphac and alphao, see
C                   subroutine AB09ID.
C
C Output parameters:
C   Ar, Br,
C   Cr, Dr - matrices of the reduced system.
C   HSV    - frequency-weighted Hankel singular values of the
C            alpha-stable part.
C   info   - warning message code:
C            info = 1 - selected order greater than the order
C                       of a minimal realization;
C            info = 2 - selected order corresponds to repeated singular
C                       values, which are neither all included nor all
C                       excluded from the reduced model;
C            info = 3 - selected order less than the order of
C                       the unstable part.
C            info = 10+K - K violations of the numerical stablity
C                       condition appeared during eigenvalue assignment.
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
      CHARACTER         DICO, EQUIL, JOB, JOBC, JOBO, ORDSEL, WEIGHT
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDV, LDW,
     $                  LDWORK, M, MW, N, NR, NS, NV, NW, P, PV
      DOUBLE PRECISION  ALPHA, ALPHAC, ALPHAO, TOL1, TOL2
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:),   B(:,:), C(:,:), D(:,:),
     $                                DWORK(:), HSV(:), V(:,:), W(:,:)
      INTEGER, ALLOCATABLE::          IWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           DISCR, LEFTW, RIGHTW
      INTEGER           I, J, LCF, LW, M1, METH, N1, N2, N3, NN, NNV,
     $                  NNW, NV1, NVP, NVP1, NW1, NWM, NWM1, P1, PPV
      DOUBLE PRECISION  ALFA(3), DUM, TOL(2)
C
C .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          LSAME, DLAMCH
C
C .. External Subroutines ..
      EXTERNAL          AB09ID, DLACPY
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FWERED requires at least 5 input arguments' )
      ELSE IF( NLHS .GT. 6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FWERED requires at most 6 output arguments' )
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
      IF( METH.LT.1 .OR. METH.GT.4 .OR. I.LT.1 .OR. I.GT.2
     $    .OR. J.LT.1 .OR. J.GT.2 ) THEN
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
C   AV(NVxNV), BV(NVxP), CV(PVxNV), DV(PVxP)
C
      NV = 0
      PV = 0
      LEFTW = .FALSE.
      IF( NRHS.GT.5 ) THEN
         NVP  = mxGetM( PRHS(6) )
         NVP1 = mxGetN( PRHS(6) )
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'V must be a real matrix' )
         END IF
         IF( NVP*NVP1.GT.0 ) THEN
            NV = NVP1 - P
            PV = NVP  - NV
            IF( NV.LT.0 .OR. PV.LT.0 ) THEN
               CALL mexErrMsgTxt( 'V incompatible with the system' )
            END IF
            LEFTW = .TRUE.
         END IF
      END IF
C
C   AW(NWxNW), BW(NWxMW), CW(MxNW), DW(MxMW)
C
      NW = 0
      MW = 0
      RIGHTW = .FALSE.
      IF( NRHS.GT.6 ) THEN
         NWM   = mxGetM( PRHS(7) )
         NWM1  = mxGetN( PRHS(7) )
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'W must be a real matrix' )
         END IF
         IF( NWM*NWM1.GT.0 ) THEN
            NW = NWM  - M
            MW = NWM1 - NW
            IF( NW.LT.0 .OR. MW.LT.0 ) THEN
               CALL mexErrMsgTxt( 'W incompatible with the system' )
            END IF
            RIGHTW = .TRUE.
         END IF
      END IF
C
C   tol(1x2)
C
      TOL1 = ZERO
      TOL2 = ZERO
      IF( NRHS.GT.7 ) THEN
         I = mxGetM( PRHS(8) )*mxGetN( PRHS(8) )
         IF( I.GT.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'TOL must be a vector with at most 2 elements' )
         END IF
         IF( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), TOL, I )
         IF( I.GT.0 ) TOL1 = TOL(1)
         IF( I.GT.1 ) TOL2 = TOL(2)
      END IF
C
C   discr
C
      DICO = 'C'
      IF( NRHS.GT.8 ) THEN
         IF( mxGetM( PRHS(9) ).NE.1 .OR. mxGetN( PRHS(9) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(9) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'DISCR must be an integer scalar 0 or 1' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), DUM, 1 )
         IF( DUM.NE.ZERO ) DICO = 'D'
      END IF
      DISCR = LSAME( DICO, 'D' )
C
C   ord
C
      ORDSEL = 'A'
      NR = 0
      IF( NRHS.GT.9 ) THEN
         IF( mxGetM( PRHS(10) ).NE.1 .OR. mxGetN( PRHS(10) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(10) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(10) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ORD must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(10) ), DUM, 1 )
         IF( DUM.GE.ZERO ) THEN
            ORDSEL = 'F'
            NR = DUM
            NR = MIN( N, NR )
         END IF
      END IF
C
C   alpha
C
      ALPHAC = ZERO
      ALPHAO = ZERO
      IF( NRHS.GT.10 ) THEN
         I = mxGetM( PRHS(11) )*mxGetN( PRHS(11) )
         IF( I.GT.3 ) THEN
            CALL mexErrMsgTxt
     $           ( 'ALPHA must be a vector with at most 3 elements' )
         END IF
         IF( mxIsNumeric( PRHS(11) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(11) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(11) ), ALFA, I )
         IF( I.GT.0 ) THEN
            ALPHA = ALFA(1)
            IF( DISCR ) THEN
               IF( ALPHA.EQ.ONE )  ALPHA = ONE - SQRT( DLAMCH( 'E' ) )
            ELSE
               IF( ALPHA.EQ.ZERO ) ALPHA = -SQRT( DLAMCH( 'E' ) )
            END IF
            IF( I.GT.1 ) ALPHAC = ALFA(2)
            IF( I.GT.2 ) ALPHAO = ALFA(3)
         ELSE
            ALPHA = -SQRT( DLAMCH( 'E' ) )
            IF( DISCR ) ALPHA = ONE + ALPHA
         END IF
      ELSE
         ALPHA = -SQRT( DLAMCH( 'E' ) )
         IF( DISCR ) ALPHA = ONE + ALPHA
      END IF
C
C Determine the lenghts of working arrays.
C
      LW  = 1
      NN  = N*N
      NNV = N + NV
      NNW = N + NW
      PPV = MAX( P, PV )
      IF( LEFTW .AND. PV.GT.0 ) THEN
         LW = MAX( LW, NNV*( NNV + MAX( NNV, PV ) + 5 ) )
      ELSE
         LW = MAX( LW, N*( P + 5 ) )
      END IF
C
      IF( RIGHTW .AND. MW.GT.0 ) THEN
         LW = MAX( LW, NNW*( NNW + MAX( NNW, MW ) + 5 ) )
      ELSE
         LW = MAX( LW, N*( M + 5 ) )
      END IF
      LW = 2*NN + MAX( LW, 2*NN + 5*N, N*MAX( M, P ) )
C
      IF( LEFTW .AND. NV.GT.0 ) THEN
         LCF = PV*( NV + PV ) + PV*NV +
     $         MAX( NV*( NV + 5 ), PV*( PV + 2 ), 4*PPV )
         IF( PV.EQ.P ) THEN
            LW = MAX( LW, LCF, NV + MAX( NV, 3*P ) )
         ELSE
            LW = MAX( LW, PPV*( 2*NV + PPV ) +
     $                    MAX( LCF, NV + MAX( NV, 3*PPV ) ) )
         END IF
      END IF
C
      IF( RIGHTW .AND. NW.GT.0 ) THEN
         IF( MW.EQ.M ) THEN
            LW = MAX( LW, NW + MAX( NW, 3*M ) )
         ELSE
            LW = MAX( LW, 2*NW*MAX( M, MW ) +
     $                    NW + MAX( NW, 3*M, 3*MW ) )
         END IF
         LW = MAX( LW, MW*( NW + MW ) +
     $             MAX( NW*( NW + 5 ), MW*( MW + 2 ), 4*MW, 4*M ) )
      END IF
      LDWORK = LW
C
      LDA = MAX( 1, N )
      LDB = MAX( 1, N )
      LDC = MAX( 1, P )
      LDD = MAX( 1, P )
      LDV = MAX( 1, NVP )
      LDW = MAX( 1, NWM )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
      ALLOCATE ( A( LDA, MAX( 1, N ) ), B( LDB, MAX( 1, M ) ),
     $           C( LDC, MAX( 1, N ) ), D( LDD, MAX( 1, M ) ),
     $           V( LDV, MAX( 1, NVP1 ) ), W( LDW, MAX( 1, NWM1 ) ),
     $           DWORK( LDWORK ), HSV( MAX( 1, N ) ),
     $           IWORK( MAX( 3, 2*N, NV + PPV, NW + MAX( MW, M ) ) ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), D, P*M )
      IF( LEFTW )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), V, NVP*NVP1 )
      IF( RIGHTW )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), W, NWM*NWM1 )
C
C Do the actual computations.
C
      EQUIL = 'S'
      IF( LEFTW .AND. RIGHTW ) THEN
         WEIGHT = 'B'
      ELSE IF( LEFTW ) THEN
         WEIGHT = 'L'
      ELSE IF( RIGHTW ) THEN
         WEIGHT = 'R'
      ELSE
         WEIGHT = 'N'
      END IF
C
      NV1 = NV + 1
      NW1 = NW + 1
      CALL AB09ID( DICO, JOBC, JOBO, JOB, WEIGHT, EQUIL, ORDSEL,
     $             N, M, P, NV, PV, NW, MW, NR, ALPHA, ALPHAC, ALPHAO,
     $             A, LDA, B, LDB, C, LDC, D, LDD,
     $             V, LDV, V(1,NV1), LDV, V(NV1,1), LDV, V(NV1,NV1),
     $             LDV, W, LDW, W(1,NW1), LDW, W(NW1,1), LDW,
     $             W(NW1,NW1), LDW, NS, HSV, TOL1, TOL2, IWORK, DWORK,
     $             LDWORK, IWARN, INFO )
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
            HSV(1)  = IWARN
            CALL mxCopyReal8ToPtr( HSV, mxGetPr( PLHS(6) ), 1 )
         END IF
C
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, D, V, W, HSV, IWORK, DWORK )
C
C Error handling.
C
      IF( INFO.NE.0 ) THEN
         WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM AB09ID" )' ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of FWERED ***
      END
