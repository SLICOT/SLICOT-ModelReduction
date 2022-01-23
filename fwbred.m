function [sysr,hsv] = fwbred(sys,varargin)
%FWBRED Frequency-weighted balancing related model reduction.
%       [SYSR,HSV] = FWBRED(SYS,V,W,TOL,ORD)  calculates for the 
%       transfer-function
%                                       -1
%               G(lambda) =  C(lambdaI-A) B + D
%
%       of an original system SYS = (A,B,C,D), an approximate 
%       transfer-function
%                                          -1
%               Gr(lambda) =  Cr(lambdaI-Ar) Br + Dr
%
%       of a reduced order system SYSR = (Ar,Br,Cr,Dr) 
%       by minimizing the frequency-weighted error norm 
%
%              ||V*(G-Gr)*W||    
%
%       using frequency-weighted balancing related approximation 
%       methods on the stable part of SYS. V and W are appropriate
%       frequency-weighting systems.
%       TOL is the tolerance for model reduction.
%       ORD specifies the desired order of the reduced system SYSR.
%
%       HSV contains the decreasingly ordered frequency-weighted Hankel 
%       singular values of the stable part of SYS.
%
%       [SYSR,HSV] = FWBRED(SYS,V,W,OPTIONS)  calculates the reduced
%       order model using the option values in the structure OPTIONS, 
%       created with the SYSREDSET function.  See SYSREDSET for details.  
%       FWBRED uses these options: BalredMethod, AccuracyEnhancing, 
%       FWEContrGramian, FWEObservGramian, TolRed, TolMinreal, 
%       CStabDeg, DStabDeg, Order.
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
%          ALPHA-stable part G1 using the frequency-weighted
%          Balance & Truncate Approximation (BTA) or Singular Perturbation
%          Approximation (SPA) method.
%   
%       3) Assemble the reduced model Gr as
%
%             Gr = G1r + G2.
%
%       RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function FWERED.
%       A. Varga 05-11-2000; revised 19-05-2001.
%       Revised, V. Sima 23-06-2001, 12-01-2002, 25-02-2009.
%

defaultopt = struct( ...
    'BalredMethod', 'bta',...
    'AccuracyEnhancing', 'bfsr', ...
    'FWEContrGramian', 'standard', ...
    'FWEObservGramian', 'standard', ...
    'TolRed', 0, ...
    'TolMinreal', 0, ...
    'CStabDeg', -sqrt(eps), ...
    'DStabDeg', 1-sqrt(eps), ...
    'FWEAlphaContr', 0, ...
    'FWEAlphaObserv', 0, ...
    'Order', -1);

% If just 'defaults' passed in, return the default options in SYSR
if nargin == 1 && nargout <= 1 && isequal(sys,'defaults')
   sysr = defaultopt;
   return
end

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end
 
discr = double(sys.ts > 0);

% initialization
ni = nargin;
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
if ni < 5 
   ord = sysredget(options,'Order',defaultopt,'fast');
else
   ord = varargin{4};
end
if ni < 4 
   tol = sysredget(options,'TolRed',defaultopt,'fast');
else
   tol = varargin{3};
end
if ni < 3 
   w = []; 
else 
   w = varargin{2};
end
if ni < 2 
   v = []; 
else
   v = varargin{1};
end

if ~isempty(v)
   if ~isa(v,'lti') && ~isa(v,'double')
      error('The left frequency weight V must be an LTI object or a matrix')
   end
   if size(v,2) ~= size(sys,1)
      error('SYS and V have incompatible dimensions')
   end
   if isa(v,'lti') && sys.ts ~= v.ts
      error('SYS and V must have the same sampling time')
   end
end
 
if ~isempty(w)
   if ~isa(w,'lti') && ~isa(w,'double')
      error('The right frequency weight W must be an LTI object or a matrix')
   end
   if size(w,1) ~= size(sys,2)
      error('SYS and W have incompatible dimensions')
   end
   if isa(w,'lti') && sys.ts ~= w.ts
      error('SYS and W must have the same sampling time')
   end
end


[a,b,c,d] = ssdata(sys);

if ~isempty(v)
   if isa(v,'double')
      vsys = v;
   else
      [av,bv,cv,dv] = ssdata(v);
      vsys = [av bv;cv dv];
   end
else
   vsys = [];
end
 
if ~isempty(w)
   if isa(w,'double')
      wsys = w;
   else
      [av,bv,cv,dv] = ssdata(w);
      wsys = [av bv;cv dv];
   end
else
   wsys = [];
end

balmeth  = sysredget(options,'BalredMethod',defaultopt,'fast');
accenh   = sysredget(options,'AccuracyEnhancing',defaultopt,'fast');
tolmin   = sysredget(options,'TolMinreal',defaultopt,'fast');
alphac   = sysredget(options,'FWEAlphaContr',defaultopt,'fast');
alphao   = sysredget(options,'FWEAlphaObserv',defaultopt,'fast');
gramc    = sysredget(options,'FWEContrGramian',defaultopt,'fast');
gramo    = sysredget(options,'FWEObservGramian',defaultopt,'fast');
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
if strcmp(gramo,'standard')
   meth = meth + 10;
else
   meth = meth + 20;
end
if strcmp(gramc,'standard')
   meth = meth + 100;
else
   meth = meth + 200;
end

[ar,br,cr,dr,hsv] = fwered(meth,a,b,c,d,vsys,wsys,[tol tolmin],discr,...
                           ord,[alpha,alphac,alphao]);

sysr = ss(ar,br,cr,dr,sys);


% end fwbred
