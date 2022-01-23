function [sysr,hsv] = spa_cf(fact,sys,tol,ord,alpha)
%SPA_CF Singular perturbation approximation of coprime factors 
%       without balancing.
%       [SYSR,HSV] = SPA_CF(FACT,SYS,TOL,ORD,ALPHA)  calculates for 
%       the transfer function
%                                    -1
%            G(lambda) =  C(lambdaI-A) B + D
%
%       of an original system SYS = (A,B,C,D) an approximate  
%       transfer function
%                                       -1
%            Gr(lambda) =  Cr(lambdaI-Ar) Br + Dr
%
%       of a reduced order system SYSR = (Ar,Br,Cr,Dr) using the
%       balancing-free Singular Perturbation Approximation method 
%       on a Left/Right Coprime Factorization (LCF/RCF) of G.
%       The allowed factorizations are:
%          FACT = 'rcfid' - RCF with inner denominator
%          FACT = 'lcfid' - LCF with inner denominator
%          FACT = 'rcfs'  - RCF with prescribed stability degree ALPHA
%          FACT = 'lcfs'  - LCF with prescribed stability degree ALPHA
%          FACT = ' '     - no coprime factorization is used
%
%       TOL is a tolerance vector [TOL1, TOL2, TOL3], where TOL1 is
%       tolerance for model reduction, TOL2 is the tolerance for
%       computing a minimal realization of the extended system
%       (see Method with 'type spabal_cf') and TOL3 is the tolerance 
%       for controlability/observability tests for computing right/left
%       coprime factorizations.
%
%       ORD specifies the desired order of the reduced system SYSR.
%
%       ALPHA is the stability degree for the eigenvalues of A when
%       using the RCF or LCF with prescribed stability degree. 
%       For a continuous-time system ALPHA <= 0 is the maximum value
%       for the real parts of eigenvalues (default: -0.05), while for
%       a discrete-time system, 1 >= ALPHA >= 0 represents the maximum
%       value for the moduli of eigenvalues (default: 0.95).
%
%       HSV contains the decreasingly ordered Hankel singular values of
%       the extended system (see Method with 'type spa_cf').
%
%       The order NR of the reduced system SYSR is determined as follows:
%       let NMIN be the order of a minimal realization of the extended
%       system (see Method with 'type spa_cf'). Then
%       (1) if TOL1 > 0 and  ORD < 0, then NR = min(NRS,NMIN), where 
%           NRS is the number of Hankel singular values greater than TOL1;
%       (2) if ORD >= 0, then NR = MIN(ORD,NMIN).
%

%       Method:
%       The following approach is used in conjunction with either
%       a stable Left Coprime Factorization (LCF) or a stable 
%       Right Coprime Factorization (RCF) of G:
%        
%       1. Compute the appropriate stable coprime factorization of G:
%                     -1                   -1
%                G = R  *Q (LCF) or G = Q*R   (RCF).
%
%       2. Perform the balancing-free Singular Perturbation Approximation
%          model reduction algorithm on the extended system
%                                           ( Q )
%                Ge = ( Q R ) (LCF) or Ge = ( R )  (RCF)
%
%          to obtain a reduced extended system with reduced factors
%                                               ( Qr )
%                Ger = ( Qr Rr ) (LCF) or Ger = ( Rr )  (RCF).
%
%       3. Recover the reduced system from the reduced factors as
%                       -1                       -1
%                Gr = Rr  *Qr (LCF) or Gr = Qr*Rr   (RCF).
%

%       RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function SYSRED.
%       A. Varga 08-03-1998. Revised 02-16-1999.
%       Revised, V. Sima 12-01-2002.
%

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end
 
ni = nargin;
discr = double(sys.ts > 0);

if ni < 5
   alpha = -0.05;
   if discr
      alpha = 1 + alpha;
   end
end
if ni < 4 
   ord = -1;
end
if ni < 3 
   tol = 0; 
end

k = find(strcmp({'rcfid','lcfid','rcfs','lcfs'},fact)==1);
if isempty(k), k = 0; end

[a,b,c,d]=ssdata(sys);

[ar,br,cr,dr,hsv]=sysred(10*k+4,a,b,c,d,tol,discr,ord,alpha);

sysr = ss(ar,br,cr,dr,sys);


% end spa_cf
