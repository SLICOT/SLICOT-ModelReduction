% SFORED.F - MEX-function for SLICOT controller reduction routines
%            SB16BD.F and SB16CD.F. 
%
% Matlab call:  
%   [Ac,Bc,Cc,Dc,HSV,info] = SFORED(meth,A,B,C,D,F,G,tol,discr,ord)
%
% Purpose:
%   To compute, for a given open-loop model (A,B,C,D), and for
%   given state feedback gain F and full observer gain G,
%   such that A+B*F and A+G*C are stable, a reduced order 
%   controller model (Ac,Bc,Cc,Dc) using a coprime factorization
%   based controller reduction approach. For reduction of
%   coprime factors, optionally a stability enforcing frequency-weighted
%   model reduction can be used.  For reduction, either the square-root
%   or the balancing-free square-root versions of the Balance & Truncate
%   (B&T) or Singular Perturbation Approximation (SPA) 
%   (only in the non-weighted case) model reduction methods are used in 
%   conjunction with stable coprime factorization techniques.
%   The order of the reduced model is determined either by the number
%   of Hankel-singular values greater than tol or by the desired
%   order ord.
%
% Description of input parameters: 
%   meth  - method flag of decimal form mfr to specify 
%           the reduction method. The allowed values for m, f and r are:
%             m = 1 : standard coprime factorization;
%             m = 2 : coprime factorization with frequency-weighting;
%             f = 1 : use left coprime factorization;
%             f = 2 : use right coprime factorization;
%             r = 1 : B&T method with balancing;
%             r = 2 : B&T method (no balancing);
%             r = 3 : SPA method with balancing (only for m = 1);
%             r = 4 : SPA (no balancing) (only for m = 1).
%   A,B,
%   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
%           and P-by-M, respectively.
%   F,G   - state-feedack and output-injection matrices of size M-by-N
%           and N-by-P, respectively.
%   tol   - (optional) tolerance vector for determining the order of 
%           reduced system, of the form [tol1, tol2], where:
%             tol1 specifies the tolerance for model reduction.
%                  Default: tol1 = N*epsilon_machine*HSV(1),
%                  where HSV(1) is the largest Hankel-singular value
%                  of the extended system Ge (see SB16BD) for m = 1 or 
%                  the largest frequecy-weighted Hankel-singular value
%                  (see SB16CD) for m = 2.
%             tol2 specifies, for m = 1, the tolerance for minimal 
%                  realization.
%                  Default: tol2 = N*epsilon_machine*HSV(1).
%   discr - (optional) type of system:
%             = 0 : continuous-time (default);
%             = 1 : discrete-time.
%   ord   - (optional) desired order of reduced system.
%             Default: ord = -1 (order determined automatically).
%
% Description of output parameters:
%   Ac, Bc, 
%   Cc, Dc - matrices of the reduced order controller.
%   HSV    - Hankel-singular values of the extended system Ge, for
%            m = 1 (see SB16BD), or the frequency-weighted 
%            Hankel-singular values, for m = 2 (see SB16CD).
%   info   - warning message code:
%            info = 1 - selected order greater than the order
%                       of a minimal realization of the controller;
%            info = 2 - selected order corresponds to repeated singular
%                       values, which are neither all included nor all 
%                       excluded from the reduced model.
%

% RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% Contributors:
%   D. Sima, University of Bucharest, and
%   A. Varga, German Aerospace Center,
%   DLR Oberpfaffenhofen, March 2001.
%
% Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, June 2001.
%
