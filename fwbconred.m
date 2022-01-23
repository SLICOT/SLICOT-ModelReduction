function [syscr,hsvc] = fwbconred(sysc,varargin)
%FWBCONRED Frequency-weighted balancing related controller reduction.
%       [SYSCR,HSVC] = FWBCONRED(SYSC,SYS,TOL,ORD)  calculates for the 
%       transfer function
%                                          -1
%               Gc(lambda) =  Cc(lambdaI-Ac) Bc + Dc
%
%       of an original controller SYSC = (Ac,Bc,Cc,Dc), an approximate 
%       transfer function
%                                             -1
%               Gcr(lambda) =  Ccr(lambdaI-Acr) Bcr + Dcr
%
%       of a reduced order controller SYSCR = (Acr,Bcr,Ccr,Dcr) by
%       minimizing the frequency-weighted error norm 
%
%              ||V*(Gc-Gcr)*W||    
%
%       using frequency-weighted balancing related approximation 
%       methods on the stable part of SYSC. V and W are special
%       frequency-weighting transfer-function matrices constructed
%       to enforce closed-loop stability and/or closed-loop performance.
%       If G is the transfer-function matrix of the open-loop system SYS, 
%       then the following weightings V and W can be used:
%                       -1
%       (a)   V = (I-G*K) *G, W = I - to enforce closed-loop stability;
%                               -1
%       (b)   V = I,  W = (I-G*K) *G - to enforce closed-loop stability;
%                       -1              -1
%       (c)   V = (I-G*K) *G, W = (I-G*K)  - to enforce closed-loop 
%             stability and performance. 
%       TOL is the tolerance for controller reduction.
%       ORD specifies the desired order of the reduced controller SYSCR.
%
%       HSVC contains the decreasingly ordered frequency-weighted Hankel 
%       singular values of the stable part of SYSC.
%
%       [SYSCR,HSVC] = FWBCONRED(SYSC,SYS,OPTIONS) calculates the 
%       reduced order controller using the option values in the structure 
%       OPTIONS, created with the SYSREDSET function.  See SYSREDSET for details.  
%       FWBCONRED uses these options: BalredMethod, AccuracyEnhancing, 
%       FWEContrGramian, FWEObservGramian, TolRed, TolMinreal, 
%       CStabDeg, DStabDeg, Order, FWEConredMethod.
%
%       The choice of frequency-weighting can be specified by the 
%       OPTIONS.FWEConredMethod structure element as follows:
%       'none'        - no weighting;
%       'outputstab'  - for the stability enforcing choice (a);
%       'inputstab'   - for the stability enforcing choice (b);
%       'performance' - for the stability and performance enforcing choice (c).
%
%       An arbitrary stability degree parameter ALPHA can be specified 
%       in the structure OPTIONS as OPTIONS.CStabDeg for a continuous-time 
%       system or OPTIONS.DStabDeg for a discrete-time system.
%       ALPHA is the stability boundary for the eigenvalues of Ac. 
%       For a continuous-time system ALPHA <= 0 is the boundary value
%       for the real parts of eigenvalues, while for a discrete-time 
%       system, 1 >= ALPHA >= 0 represents the boundary value for the 
%       moduli of eigenvalues.
%
%       The order NR of the reduced controller SYSCR is determined as follows:
%       let NU be the order of the ALPHA-unstable part of SYSC and let 
%       NSMIN be the order of a minimal realization of the ALPHA-stable 
%       part. Then
%       (1) if TOL > 0 and  ORD < 0, then NR = NU + min(NRS,NSMIN), where 
%           NRS is the number of Hankel singular values greater than TOL;
%       (2) if ORD >= 0, then NR = NU + min(max(0,ORD-NU),NSMIN).
%

%       Method: 
%       The following approach is used to reduce a given Gc:
%
%       1) Decompose additively Gc as
%
%            Gc = Gc1 + Gc2
%
%          such that Gc1 = (Acs,Bcs,Ccs,Dc) has only ALPHA-stable poles and 
%          Gc2 = (Acu,Bcu,Ccu,0) has only ALPHA-unstable poles.
%
%       2) Determine Gc1r, a reduced order approximation of the 
%          ALPHA-stable part Gc1 using the frequency-weighted
%          Balance & Truncate Approximation method.
%   
%       3) Assemble the reduced model Gcr as
%
%             Gcr = Gc1r + Gc2.
%
%       RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function CONRED.
%       A. Varga 05-11-2000; revised 22-05-2001.
%       Revised, V. Sima 23-06-2001, 12-01-2002, 25-02-2009.
%

defaultopt = struct( ...
    'BalredMethod', 'bta',...
    'AccuracyEnhancing', 'bfsr', ...
    'FWEContrGramian', 'standard', ...
    'FWEObservGramian', 'standard', ...
    'FWEConredMethod', 'performance', ...
    'TolRed', 0, ...
    'TolMinreal', 0, ...
    'CStabDeg', -sqrt(eps), ...
    'DStabDeg', 1-sqrt(eps), ...
    'FWEAlphaContr', 0, ...
    'FWEAlphaObserv', 0, ...
    'Order', -1);
 
% If just 'defaults' passed in, return the default options in SYSCR
if nargin == 1 && nargout <= 1 && isequal(sysc,'defaults')
   syscr = defaultopt;
   return
end

if ~isa(sysc,'lti')
   error('The input controller SYSC must be an LTI object')
end
 
 
ni = nargin;
discr = double(sysc.ts > 0);

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
if ni < 4 
   ord = sysredget(options,'Order',defaultopt,'fast');
else
   ord = varargin{3};
end
if ni < 3 
   tol = sysredget(options,'TolRed',defaultopt,'fast');
else
   tol = varargin{2};
end
if ni < 2 
   sys = []; 
else
   sys = varargin{1};
end

if ~isempty(sys)
   if ~isa(sys,'lti') 
      error('SYS must be an LTI object or a matrix')
   end
   if size(sys,2) ~= size(sysc,1) || size(sys,1) ~= size(sysc,2)
      error('SYSC and SYS have incompatible dimensions')
   end
   if sysc.ts ~= sys.ts
      error('SYSC and SYS must have the same sampling time')
   end
   [a,b,c,d] = ssdata(sys);
else
    a = []; b = []; c = []; d = [];
end
 


[ac,bc,cc,dc] = ssdata(sysc);

balmeth    = sysredget(options,'BalredMethod',defaultopt,'fast');
accenh     = sysredget(options,'AccuracyEnhancing',defaultopt,'fast');
tolmin     = sysredget(options,'TolMinreal',defaultopt,'fast');
gramc      = sysredget(options,'FWEContrGramian',defaultopt,'fast');
gramo      = sysredget(options,'FWEObservGramian',defaultopt,'fast');
conredmeth = sysredget(options,'FWEConredMethod',defaultopt,'fast');
if discr
   alpha   = sysredget(options,'DStabDeg',defaultopt,'fast');
else
   alpha   = sysredget(options,'CStabDeg',defaultopt,'fast');
end

if strcmp(conredmeth,'none')
   meth = 1;
elseif strcmp(conredmeth,'outputstab')
   meth = 2;
elseif strcmp(conredmeth,'inputstab')
   meth = 3;
else
   meth = 4;
end
if strcmp(balmeth,'bta')
   meth = meth + 10;
else
   meth = meth + 30;
end
if strcmp(accenh,'bfsr')
   meth = meth + 10;
end
if strcmp(gramo,'standard')
   meth = meth + 100;
else
   meth = meth + 200;
end
if strcmp(gramc,'standard')
   meth = meth + 1000;
else
   meth = meth + 2000;
end

[acr,bcr,ccr,dcr,hsvc] = conred(meth,ac,bc,cc,dc,a,b,c,d,[tol tolmin],discr,...
                                ord,alpha);

syscr = ss(acr,bcr,ccr,dcr,sysc);


% end fwbconred
