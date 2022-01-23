C FWEHNA.F - Gateway function for SLICOT model reduction routine
C            AB09JD.F.
C
C RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
C Based on SLICOT RELEASE 5.7. Copyright (c) 2002-2020 NICONET e.V.
C
C Matlab call:
C   [Ar,Br,Cr,Dr,HSV,info] = FWEHNA(meth,A,B,C,D,V,W,tol,discr,ord,alpha)
C
C Purpose:
C   To compute a reduced order model (Ar,Br,Cr,Dr) for an original
C   state-space representation (A,B,C,D) by using the frequency
C   weighted optimal Hankel-norm approximation method.
C   The Hankel norm of the weighted error
C
C         op(V)*(G-Gr)*op(W)
C
C   is minimized, where G and Gr are the transfer-function matrices
C   of the original and reduced systems, respectively, V and W
C   are invertible transfer-function matrices of the left and right
C   frequency weights, and op(X) denotes X, conj(X), inv(X), or
C   conj(inv(X)). V and W are specified by their state space
C   realizations (AV,BV,CV,DV) and (AW,BW,CW,DW), respectively.
C   When minimizing the weighted error V*(G-Gr)*W, V and W must have
C   poles distinct from those of G.
C   When minimizing conj(V)*(G-Gr)*conj(W), conj(V) and conj(W) must
C   have poles distinct from those of G.
C   Additionally, V and W must be invertible transfer-function matrices.
C   If the original system is unstable, then the frequency-weighted
C   Hankel-norm approximation is computed only for the
C   alpha-stable part of the system.
C   The order of the reduced model is determined either by the number
C   of frequency-weighted Hankel-singular values HSV greater than tol or
C   by the desired order ord.
C
C Input parameters:
C   meth  - method flag of decimal form ijk, where:
C             i = 1 : op(V) = V;
C             i = 2 : op(V) = conj(V);
C             i = 3 : op(V) = inv(V);
C             i = 4 : op(V) = conj(inv(V));
C             j = 1 : op(W) = W;
C             j = 2 : op(W) = conj(W);
C             j = 3 : op(W) = inv(W);
C             j = 4 : op(W) = conj(inv(W));
C             k = 1 : if possible, use standard inverse based method;
C             k = 2 : use inverse free method;
C             k = 3 : use an automatic method.
C   A,B,
C   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
C           and P-by-M, respectively.
C   V     - (optional) contains the system matrix [AV BV; CV DV] of
C           the left weighting.
C   W     - (optional) contains the system matrix [AW BW; CW DW] of
C           the right weighting.
C   tol   - (optional) tolerance vector for determining the order of
C           reduced system, of the form [tol1, tol2], where:
C             tol1 specifies the tolerance for model reduction.
C                  Default: tol1 = NS*epsilon_machine*Hankel_norm(G1s),
C                  where Hankel_norm(G1s) is the Hankel-norm of the
C                  projection G1s of op(V)*G1*op(W) (see AB09JD,
C                  Section METHOD) and NS is the order of G1s.
C             tol2 specifies the tolerance for computing a minimal
C                  realization of the alpha-stable part of the
C                  weighted original system.
C                  Default: tol2 = NS*epsilon_machine*Hankel_norm(G1s).
C   discr - (optional) type of system:
C              = 0 : continuous-time (default);
C              = 1 : discrete-time.
C   ord   - (optional) desired order of reduced system.
C             Default: ord = -1 (order determined automatically).
C   alpha - (optional) stability boundary for the eigenvalues of A.
C             Default:    -sqrt(epsilon_machine)  for continuous-time;
C                      1.0-sqrt(epsilon_machine)  for discrete-time.
C
C Output parameters:
C   Ar, Br,
C   Cr, Dr - state-space matrices of the reduced system Gr.
C   HSV    - Hankel singular values of the projection G1s of
C            op(V)*G1*op(W) (see AB09JD, Section METHOD),
C            where G1 is the alpha-stable part of G.
C   info   - warning message code:
C            info = 1 - selected order greater than the order
C                       of a minimal realization;
C            info = 2 - selected order less than the order of
C                       the unstable part.
C
C Contributors:
C   D. Sima, University of Bucharest, and
C   A. Varga, German Aerospace Center,
C   DLR Oberpfaffenhofen, March 2001.
C
C Revisions:
C   V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C   March 2005, Apr. 2009, Dec. 2012.
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
      CHARACTER         DICO, EQUIL, JOBINV, JOBV, JOBW, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDV, LDW,
     $                  LDWORK, M, N, NR, NS, NV, NW, P
      DOUBLE PRECISION  ALPHA, TOL1, TOL2
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
      INTEGER           I, J, LW, LWI, M1, METH, N1, N2, N3,
     $                  NV1, NVP, NVP1, NW1, NWM, NWM1, P1
      DOUBLE PRECISION  DUM, TOL(2)
C
C .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C
C .. External Subroutines ..
      EXTERNAL          AB09JD, DLACPY
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FWEHNA requires at least 5 input arguments' )
      ELSE IF( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'FWEHNA requires at most 6 output arguments' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   meth = i*100+j*10+k
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
      IF( METH.LT.1 .OR. METH.GT.3 .OR. I.LT.1 .OR. I.GT.4
     $    .OR. J.LT.1 .OR. J.GT.4 ) THEN
         CALL mexErrMsgTxt( 'Invalid value for METH ' )
      END IF
      IF( I.EQ.1 ) THEN
         JOBV = 'V'
      ELSE IF( I.EQ.2 ) THEN
         JOBV = 'C'
      ELSE IF( I.EQ.3 ) THEN
         JOBV = 'I'
      ELSE
         JOBV = 'R'
      END IF
      IF( J.EQ.1 ) THEN
         JOBW = 'W'
      ELSE IF( J.EQ.2 ) THEN
         JOBW = 'C'
      ELSE IF( J.EQ.3 ) THEN
         JOBW = 'I'
      ELSE
         JOBW = 'R'
      END IF
      IF( METH.EQ.1 ) THEN
         JOBINV = 'I'
      ELSE IF( METH.EQ.2 ) THEN
         JOBINV = 'N'
      ELSE
         JOBINV = 'A'
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
C   AV(NVxNV), BV(NVxP), CV(PxNV), DV(PxP)
C
      NV = 0
      LEFTW = .FALSE.
      IF( NRHS.GT.5 ) THEN
         NVP  = mxGetM( PRHS(6) )
         NVP1 = mxGetN( PRHS(6) )
         IF( NVP1.NE.NVP ) THEN
            CALL mexErrMsgTxt
     $          ( 'V must be a square matrix' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'V must be a real matrix' )
         END IF
         IF( NVP.GT.0 ) THEN
            NV = NVP - P
            IF( NV.LT.0 ) THEN
               CALL mexErrMsgTxt( 'V incompatible with the system' )
            END IF
         END IF
         LEFTW = NVP.GT.0
      END IF
      IF( .NOT.LEFTW ) JOBV = 'N'
C
C   AW(NWxNW), BW(NWxM), CW(MxNW), DW(MxM)
C
      NW = 0
      RIGHTW = .FALSE.
      IF( NRHS.GT.6 ) THEN
         NWM  = mxGetM( PRHS(7) )
         NWM1 = mxGetN( PRHS(7) )
         IF( NWM1.NE.NWM ) THEN
            CALL mexErrMsgTxt
     $          ( 'W must be a square matrix' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'W must be a real matrix' )
         END IF
         IF( NWM.GT.0 ) THEN
            NW = NWM - M
            IF( NW.LT.0 ) THEN
               CALL mexErrMsgTxt( 'W incompatible with the system' )
            END IF
         END IF
         RIGHTW = NWM.GT.0
      END IF
      IF( .NOT.RIGHTW ) JOBW = 'N'
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
      IF( NRHS.GT.10 ) THEN
         IF( mxGetM( PRHS(11) ).NE.1 .OR. mxGetN( PRHS(11) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(11) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(11) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(11) ), ALPHA, 1 )
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
      LWI = MAX( 1, M )
      IF( DISCR )
     $   LWI = MAX( LWI, N )
      IF( LEFTW ) THEN
         NVP = NV + P
         LW  = MAX( LW, 2*NVP*( NVP + P ) + P*P +
     $              MAX( 2*NVP*NVP + MAX( 11*NVP + 16, P*NVP ),
     $                   NVP*N + MAX( NVP*N+N*N, P*N, P*M ) ) )
         LWI = MAX( LWI, 2*P, NVP+N+6, 2*NV+P+2 )
      END IF
      IF( RIGHTW ) THEN
         NWM = NW + M
         LW  = MAX( LW, 2*NWM*( NWM + M ) + M*M +
     $              MAX( 2*NWM*NWM + MAX( 11*NWM + 16, M*NWM ),
     $                   NWM*N + MAX( NWM*N+N*N, M*N, P*M ) ) )
         LWI = MAX( LWI, 2*M, NWM+N+6, 2*NW+M+2 )
      END IF
      LW = MAX( LW, N*( 2*N + MAX( N, M, P ) + 5 ) + ( N*( N + 1 ) )/2 )
      LW = MAX( LW, N*( M + P + 2 ) + 2*M*P + MIN( N, M ) +
     $                             MAX ( 3*M + 1, MIN( N, M ) + P ) )
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
     $           V( LDV, MAX( 1, NVP ) ), W( LDW, MAX( 1, NWM ) ),
     $           DWORK( LDWORK ), HSV( MAX( 1, N ) ),
     $           IWORK( LWI ) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, N*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), D, P*M )
      IF( NVP.GT.0 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), V, NVP*NVP )
      IF( NWM.GT.0 )
     $   CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), W, NWM*NWM )
C
C Do the actual computations.
C
      EQUIL = 'S'
C
C Frequency-weighted Hankel-Norm Approximation.
C
      NV1 = NV + 1
      NW1 = NW + 1
      CALL AB09JD( JOBV, JOBW, JOBINV, DICO, EQUIL, ORDSEL, N, NV, NW,
     $             M, P, NR, ALPHA, A, LDA, B, LDB, C, LDC, D, LDD,
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
            HSV(1) = IWARN
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
         WRITE( TEXT,'( " INFO = ", I4, " ON EXIT FROM AB09JD" )' ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of FWEHNA ***
      END
