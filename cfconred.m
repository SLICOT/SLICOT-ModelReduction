function [syscr,hsv] = cfconred(f,g,sys,varargin)
%CFCONRED Coprime factorization based reduction of state feedback controllers.
%       [SYSCR,HSV] = CFCONRED(F,G,SYS,TOL,ORD)  calculates for the 
%       transfer function
%                                          -1
%               Gc(lambda) =  Cc(lambdaI-Ac) Bc + Dc
%
%       of an original state-feedback controller 
%       SYSC = (A+B*F+C*G+F*D*G,G,F,0), an approximate 
%       transfer function
%                                             -1
%               Gcr(lambda) =  Ccr(lambdaI-Acr) Bcr + Dcr
%
%       of a reduced order controller SYSCR = (Acr,Bcr,Ccr,Dcr) using
%       coprime factorization techniques in conjunction with 
%       balancing related approximation methods. The state feedback gain F and
%       the full observer gain G must be such that A+B*F and A+G*C are stable.
%       TOL is the tolerance for controller reduction.
%       ORD specifies the desired order of the reduced controller SYSCR.
%
%       HSV contains the decreasingly ordered Hankel singular values
%       defined according to the used reduction approach.
%
%       [SYSCR,HSV] = CFCONRED(F,G,SYS,OPTIONS) calculates the 
%       reduced order controller using the option values in the structure 
%       OPTIONS, created with the SYSREDSET function.  See SYSREDSET for details.  
%       CFCONRED uses these options: BalredMethod, AccuracyEnhancing, 
%       TolRed, TolMinreal, Order, CFConredMethod, CoprimeFactorization.
%
%       The choice of reduction approach can be specified by the 
%       OPTIONS.CFConredMethod structure element as follows:
%       'fwe'   - frequency-weighting (default);
%       'nofwe' - no frequency-weighting.
%    
%       The type of coprime factorization can be specified by the 
%       OPTIONS.CoprimeFactorization structure element as follows:
%       'left'  - left coprime factorization;
%       'right' - right coprime factorization (default).
%
%       The order NR of the reduced controler SYSCR is determined as follows:
%       let NSMIN be the order of the controller if frequency-weighting is used
%       or the order of a minimal realization of the extended system
%       constructed from the coprime factors of the controller. Then
%       (1) if TOL > 0 and  ORD < 0, then NR = min(NRS,NSMIN), where 
%           NRS is the number of Hankel singular values greater than TOL;
%       (2) if ORD >= 0, then NR = min(ORD,NSMIN).
%

%       RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%       Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%       Interface M-function to the SLICOT-based MEX-function SFORED.
%       A. Varga 05-11-2000; revised 24-05-2001.
%       Revised, V. Sima 23-06-2001, 12-01-2002, 25-02-2009.
%

defaultopt = struct( ...
    'BalredMethod', 'bta',...
    'AccuracyEnhancing', 'bfsr', ...
    'CoprimeFactorization', 'right', ...
    'CFConredMethod', 'fwe', ...
    'TolRed', 0, ...
    'TolMinreal', 0, ...
    'Order', -1);
 
% If just 'defaults' passed in, return the default options in SYSCR
if nargin == 1 && nargout <= 1 && isequal(f,'defaults')
   syscr = defaultopt;
   return
end

% initialization
ni = nargin;
if ni < 3
   error(' At least three arguments must be provided')
end

if nargin > 3
   if isstruct(varargin{nargin-3})
      options = varargin{nargin-3};
      ni = ni - 1;
   else
      options = []; 
   end
else
   options = []; 
end
if ni < 5 
   ord = sysredget(options,'Order',defaultopt,'fast');
else
   ord = varargin{2};
end
if ni < 4 
   tol = sysredget(options,'TolRed',defaultopt,'fast');
else
   tol = varargin{1};
end

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI object')
end
discr = double(sys.ts > 0);
[a,b,c,d] = ssdata(sys);


balmeth = sysredget(options,'BalredMethod',defaultopt,'fast');
accenh  = sysredget(options,'AccuracyEnhancing',defaultopt,'fast');
tolmin  = sysredget(options,'TolMinreal',defaultopt,'fast');
cof     = sysredget(options,'CoprimeFactorization',defaultopt,'fast');

cfconredmeth = sysredget(options,'CFConredMethod',defaultopt,'fast');

if strcmp(balmeth,'bta')
   meth = 1;
else
   meth = 3;
end
if strcmp(accenh,'bfsr')
   meth = meth + 1;
end
if strcmp(cof,'left')
   meth = meth + 10;
else
   meth = meth + 20;
end
if strcmp(cfconredmeth,'nofwe')
   meth = meth + 100;
else
   meth = meth + 200;
end

[acr,bcr,ccr,dcr,hsv] = sfored(meth,a,b,c,d,f,g,[tol tolmin],discr,...
                               ord);

syscr = ss(acr,bcr,ccr,dcr,sys);


% end cfconred
