function [sysr,hsv] = bst(sys,varargin)
%BST    Balanced stochastic truncation approximation.
%       [SYSR,HSV] = BST(SYS,TOL,ORD)  calculates for the 
%       transfer function matrix
%                                       -1
%               G(lambda) =  C(lambdaI-A) B + D
%
%       of an original system SYS = (A,B,C,D), an approximate 
%       transfer function matrix
%                                          -1
%               Gr(lambda) =  Cr(lambdaI-Ar) Br + Dr
%
%       of a reduced order system SYSR = (Ar,Br,Cr,Dr) using the
%       Balanced Stochastic Truncation (BST) approximation method on 
%       the stable part of SYS. 
%       TOL is the tolerance for model reduction.
%       ORD specifies the desired order of the reduced system SYSR.
%
%       HSV contains the decreasingly ordered frequency-weighted Hankel 
%       singular values of the stable part of SYS.
%
%       [SYSR,HSV] = BST(SYS,OPTIONS)  calculates the reduced 
%       order model using the option values in the structure OPTIONS, 
%       created with the SYSREDSET function.  See SYSREDSET for details.  
%       BST uses these options: BalredMethod, AccuracyEnhancing, 
%       TolRed, TolMinreal, CStabDeg, DStabDeg, BstBeta, Order.
%
%       An arbitrary stability degree parameter ALPHA can be specified 
%       in the structure OPTIONS as OPTIONS.CStabDeg for a continuous-time 
%       system or OPTIONS.DStabDeg for a discrete-time system.
%       ALPHA is the stability boundary for the eigenvalues of A. 
%       For a continuous-time system ALPHA <= 0 is the boundary value
%       for the real parts of eigenvalues, while for a discrete-time 
%       system, 1 >= ALPHA >= 0 represents the boundary value for the 
%       moduli of eigenvalues.
%
%       A relative-absolute weighting parameter BETA can be specified in the
%       structure OPTIONS as OPTIONS.BSTBeta. If BETA > 0, the reduction
%       is performed on the extended system [G  BETA*I]. For BETA = 0,
%       a pure relative method is employed (BST), but for BETA > 0,
%       a weighted method results. For a very large BETA the BST method
%       is practically equivalent to the Balanced Truncation method
%       applied to G. 
%       NOTE: If BETA = 0, D must have full row or column rank. If D has not
%             a full rank, then BETA = 0.01 is employed.
%
%       The order NR of the reduced system SYSR is determined as follows:
%       let NU be the order of the ALPHA-unstable part of SYS and let 
%       NSMIN be the order of a minimal realization of the ALPHA-stable 
%       part. Then
%       (1) if TOL > 0 and  ORD < 0, then NR = NU + min(NRS,NSMIN), where 
%           NRS is the number of Hankel singular values greater than TOL;
%       (2) if ORD >= 0, then NR = NU + min(max(0,ORD-NU),NSMIN).
%

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
%          stable part G1 using the BST approximation method.
%   
%       3) Assemble the reduced model Gr as
%
%             Gr = G1r + G2.
%

%       RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function BSTRED.
%       A. Varga 27-11-2000; revised 19-05-2001. 
%       Revised, V. Sima 23-06-2001, 12-01-2002, 25-02-2009.
%

defaultopt = struct( ...
    'BalredMethod','bta',...
    'AccuracyEnhancing', 'bfsr', ...
    'TolRed', 0, ...
    'TolMinreal', 0, ...
    'CStabDeg', -sqrt(eps), ...
    'DStabDeg', 1-sqrt(eps), ...
    'BstBeta',0, ...
    'Order', -1);
 
% If just 'defaults' passed in, return the default options in SYSR
if nargin == 1 && nargout <= 1 && isequal(sys,'defaults')
   sysr = defaultopt;
   return
end

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end
 
ni = nargin;
discr = double(sys.ts > 0);

% initialization
if nargin > 1
   if isstruct(varargin{nargin-1})
      options = varargin{nargin-1};
      ni = ni-1;
   else
      options = []; 
   end
else
   options = []; 
end
if ni < 3 
   ord = sysredget(options,'Order',defaultopt,'fast');
else
   ord = varargin{2};
end
if ni < 2 
   tol = sysredget(options,'TolRed',defaultopt,'fast');
else
   tol = varargin{1};
end

dualsys = 0;
beta = sysredget(options,'BstBeta',defaultopt,'fast');
if beta == 0 
   %  arrange that G has full row rank, as needed by 'bstred'
   [p,m] = size(sys.d);
   if p > m
      sys = sys.'; dualsys = 1; p = m; 
   end
   if rank(sys.d) < p, 
      beta = 0.01; 
   end
end
    
balmeth  = sysredget(options,'BalredMethod',defaultopt,'fast');
accenh   = sysredget(options,'AccuracyEnhancing',defaultopt,'fast');
tolmin   = sysredget(options,'TolMinreal',defaultopt,'fast');
if discr
   alpha = sysredget(options,'DStabDeg',defaultopt,'fast');
else
   alpha = sysredget(options,'CStabDeg',defaultopt,'fast');
end

if strcmp(balmeth,'bta')
   meth = 1;
else
   meth = 3;
end
if strcmp(accenh,'bfsr')
   meth = meth + 1;
end

[a,b,c,d] = ssdata(sys);

[ar,br,cr,dr,hsv] = bstred(meth,a,b,c,d,[tol,tolmin],discr,ord,alpha,beta);

sysr = ss(ar,br,cr,dr,sys);

if dualsys 
   sysr = sysr.'; 
end


% end bst
