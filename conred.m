% CONRED.F - MEX-function for SLICOT controller reduction routine
%            SB16AD.F. 
%
% Matlab call:  
%   [Acr,Bcr,Ccr,Dcr,HSVC,info] = CONRED(meth,Ac,Bc,Cc,Dc,A,B,C,D,...
%                                        tol,discr,ord,alpha)
%
% Purpose:
%   To compute a reduced order controller (Acr,Bcr,Ccr,Dcr) for an 
%   original state-space controller representation (Ac,Bc,Cc,Dc) by 
%   using the frequency-weighted square-root or balancing-free 
%   square-root Balance & Truncate (B&T) or Singular Perturbation 
%   Approximation (SPA) model reduction methods. The algorithm tries 
%   to minimize the norm of the frequency-weighted error
%
%           ||V*(K-Kr)*W||    
%     
%   where K and Kr are the transfer-function matrices of the original
%   and reduced order controllers, respectively. V and W are special
%   frequency-weighting transfer-function matrices constructed
%   to enforce closed-loop stability and/or closed-loop performance.
%   If G is the transfer-function matrix of the open-loop system, then
%   the following weightings V and W can be used:
%                      -1
%      (a)   V = (I-G*K) *G, W = I - to enforce closed-loop stability;
%                              -1
%      (b)   V = I,  W = (I-G*K) *G - to enforce closed-loop stability;
%                      -1              -1
%      (c)   V = (I-G*K) *G, W = (I-G*K)  - to enforce closed-loop 
%            stability and performance. 
%
%   G has the state space representation (A,B,C,D).   
%   If K is unstable, only the ALPHA-stable part of K is reduced.
%
% Description of input parameters: 
%   meth  - method flag of decimal form ijkl, where
%             i = 1 : use standard choice for controllability Grammian;
%             i = 2 : use stability garanteeing choice for the 
%                     controllability Grammian;
%             j = 1 : use standard choice for observability Grammian;
%             j = 2 : use stability garanteeing choice for the 
%                     observability Grammian;
%             k = 1 : use the square-root BT method;
%             k = 2 : use the balancing-free square-root BT method;
%             k = 3 : use the square-root SPA method;
%             k = 4 : use the balancing-free square-root SPA method;
%             l = 1 : no weightings are used;
%             l = 2 : stability enforcing left (output) weighting;
%             l = 3 : stability enforcing right (input) weighting;
%             l = 4 : stability and performance enforcing weightings.
%           Note: For a complete explanation on Grammian choices 
%                 see subroutine SB16AD.
%   AC,BC,
%   CC,DC - state-space system matrices of size NC-by-NC, NC-by-P, 
%           M-by-NC, and M-by-P, respectively.
%   A,B,
%   C,D   - (optional) state-space system matrices of size 
%           N-by-N, N-by-M, P-by-N, and P-by-M, respectively.
%   tol   - (optional) tolerance vector for determining the order of 
%           reduced system, of the form [tol1, tol2], where:
%             tol1 specifies the tolerance for model reduction.
%                  Default: tol1 = NCS*epsilon_machine*HSVC(1), where
%                  NCS is the order of the alpha-stable part of K.
%             tol2 specifies the tolerance for minimal realization.
%                  Default: tol2 = NCS*epsilon_machine*HSVC(1).
%   discr - (optional) type of system:
%              = 0 : continuous-time (default);
%              = 1 : discrete-time.
%   ord   - (optional) desired order of reduced system.
%              Default: ord = -1 (order determined automatically).
%   alpha - (optional) stability boundary for the eigenvalues of AC.
%              Default:    -sqrt(epsilon_machine)  for continuous-time;
%                       1.0-sqrt(epsilon_machine)  for discrete-time.
%
% Description of output parameters:
%   Acr, Bcr, 
%   Ccr, Dcr - matrices of the reduced controller.
%   HSVC   - frequency-weighted Hankel singular values of the 
%            alpha-stable part of K.
%   info   - warning message code:
%            info = 1 - selected order greater than the order
%                       of a minimal realization;
%            info = 2 - selected order corresponds to repeated singular
%                       values, which are neither all included nor all 
%                       excluded from the reduced model;
%            info = 3 - selected order less than the order of
%                       the unstable part.
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
