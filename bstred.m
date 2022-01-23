% BSTRED.F - MEX-function for SLICOT model reduction routine
%            AB09HD.F.
%
% Matlab call:  
%   [Ar,Br,Cr,Dr,HSV,info] = BSTRED(meth,A,B,C,D,tol,discr,ord,alpha,beta)
%
% Purpose:
%   To find a reduced order state-space system Gr = (Ar,Br,Cr,Dr)
%   from a continuous- or discrete-time original system G = (A,B,C,D) 
%   using the balanced stochastic truncation (BST) or the balanced 
%   stochastic singular perturbation approximation (BS-SPA) methods.
%   The order of the reduced model is determined either by the number
%   of stochastic Hankel-singular values HSV greater than tol or 
%   by the desired order ord.
%
% Description of input parameters: 
%   meth  - method flag to specify the basic model reduction method;
%           Allowed values for meth are:
%             meth = 1 : BST method with balancing;
%             meth = 2 : BST method (no balancing);
%             meth = 3 : BS-SPA method with balancing;
%             meth = 4 : BS-SPA (no balancing).
%   A,B,
%   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
%           and P-by-M, respectively.
%   tol   - (optional) tolerance vector for determining the order of 
%           reduced system, of the form [tol1, tol2], where:
%             tol1 specifies the tolerance for model reduction.
%                  Default: tol1 = NS*epsilon_machine, where NS is the
%                  order of the alpha-stable part of G.
%             tol2 specifies the tolerance for computing a minimal 
%                  realization when meth = 3 or 4.
%                  Default: tol2 = NS*epsilon_machine.
%   discr - (optional) type of system:
%              = 0 : continuous-time (default);
%              = 1 : discrete-time.
%   ord   - (optional) desired order of reduced system.
%             Default: ord = -1 (order determined automatically).
%   alpha - (optional) stability boundary for the eigenvalues of A.
%             Default:    -sqrt(epsilon_machine)  for continuous-time;
%                      1.0-sqrt(epsilon_machine)  for discrete-time.
%   beta  - (optional) absolute/relative error weighting parameter.
%           beta must be positive if D has not a full row rank.
%             Default: 0 (pure relative method).
%
% Description of output parameters:
%   Ar, Br, 
%   Cr, Dr - matrices of the reduced system.
%   HSV    - Hankel singular values of the alpha-stable part.
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
