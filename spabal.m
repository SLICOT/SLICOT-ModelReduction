function [sysr,hsv] = spabal(sys,tol,ord,alpha)
%SPABAL Singular perturbation approximation with balancing.
%       [SYSR,HSV] = SPABAL(SYS,TOL,ORD,ALPHA)  calculates for the 
%       transfer function
%                                       -1
%               G(lambda) =  C(lambdaI-A) B + D
%
%       of an original system SYS = (A,B,C,D) an approximate 
%       transfer function
%                                          -1
%               Gr(lambda) =  Cr(lambdaI-Ar) Br + Dr
%
%       of a reduced order system SYSR = (Ar,Br,Cr,Dr) using the
%       Singular Perturbation Approximation method on the ALPHA-stable
%       part of SYS (see Method with 'type btabal'). For a stable 
%       system SYS, the resulting reduced system SYSR is balanced.
%
%       TOL is a tolerance vector [TOL1, TOL2], where TOL1 is the  
%       tolerance for model reduction and TOL2 is the tolerance 
%       for computing a minimal realization of the ALPHA-stable part.
%
%       ORD specifies the desired order of the reduced system SYSR.
%
%       ALPHA is the stability boundary for the eigenvalues of A. 
%       For a continuous-time system ALPHA <= 0 is the boundary value
%       for the real parts of eigenvalues, while for a discrete-time 
%       system, 1 >= ALPHA >= 0 represents the boundary value for the 
%       moduli of eigenvalues.
%
%       HSV contains the decreasingly ordered Hankel singular values of
%       the ALPHA-stable part of SYS.
%
%       The order NR of the reduced system SYSR is determined as follows:
%       let NU be the order of the ALPHA-unstable part of SYS and let 
%       NSMIN be the order of a minimal realization of the ALPHA-stable 
%       part (number of Hankel singular values greater than TOL2). Then
%       (1) if TOL1 > 0 and  ORD < 0, then NR = NU + min(NRS,NSMIN), where 
%           NRS is the number of Hankel singular values greater than TOL1;
%       (2) if ORD >= 0, then NR = NU+MIN(MAX(0,ORD-NU),NSMIN).

%       Method: 
%       The following approach is used to reduce a given G:
%
%       1) Decompose additively G as
%
%            G = G1 + G2
%
%          such that G1 = (As,Bs,Cs,D) has only ALPHA-stable poles and 
%          G2 = (Au,Bu,Cu,0) has only ALPHA-unstable poles.
%
%       2) Determine G1r, a reduced order approximation of the 
%          ALPHA-stable part G1 using the Singular Perturbation 
%          Approximation method with balancing.
%   
%       3) Assemble the reduced model Gr as
%
%             Gr = G1r + G2.
%
%       RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function SYSRED.
%       A. Varga 04-05-1998. Revised 04-12-1998, 02-16-1999.
%       Revised, V. Sima 12-01-2002.
%

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end
 
ni = nargin;
discr = double(sys.ts > 0);

if ni < 4
   alpha = -sqrt(eps);
   if discr
      alpha = 1 + alpha;
   end
end
if ni < 3 
   ord = -1;
end
if ni < 2 
   tol = 0; 
end

[a,b,c,d]=ssdata(sys);

[ar,br,cr,dr,hsv]=sysred(3,a,b,c,d,tol,discr,ord,alpha);

sysr = ss(ar,br,cr,dr,sys);


% end spabal
