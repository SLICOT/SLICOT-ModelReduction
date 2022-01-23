% FWERED.F - MEX-function for SLICOT model reduction routine
%            AB09ID.F. 
%
% Matlab call:  
%   [Ar,Br,Cr,Dr,HSV,info] = FWERED(meth,A,B,C,D,V,W,tol,discr,ord,alpha)
%
% Purpose:
%   To find for a given state-space system (A,B,C,D) with 
%   transfer-function matrix G a reduced order state-space 
%   representation (Ar,Br,Cr,Dr) with the transfer-function matrix Gr
%   such that the frequency-weighted approximation error V*(G-Gr)*W is
%   minimized, where V and W are stable frequency-weighting 
%   transfer-function matrices. 
%   The frequency-weighted Hankel singular values HSV are also computed.
%   If G is unstable, the approximation is only performed for the  
%   alpha-stable part of G. 
%   The order of the reduced model is determined either by the number
%   of Hankel-singular values greater than tol or by the desired
%   order ord.
%   A frequency-weighted Balance & Truncate (BT) or 
%   Singular Perturbation Approximation (SPA) approach is used.
%
% Description of input parameters: 
%   meth  - method flag of decimal form ijk, where
%             i = 1 : use standard choice for controllability Grammian;
%             i = 2 : use stability garanteeing choice for the 
%                     controllability Grammian;
%             j = 1 : use standard choice for observability Grammian;
%             j = 2 : use stability garanteeing choice for the 
%                     observability Grammian;
%             k = 1 : use the square-root BT method;
%             k = 2 : use the balancing-free square-root BT method;
%             k = 3 : use the square-root SPA method;
%             k = 4 : use the balancing-free square-root BT method.
%           Note: For a complete explanation on Grammian choices 
%                 see subroutine AB09ID.
%   A,B,
%   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
%           and P-by-M, respectively.
%   V     - (optional) contains the system matrix [AV BV; CV DV] of 
%           the left weighting; default: V = I.
%   W     - (optional) contains the system matrix [AW BW; CW DW] of 
%           the right weighting; default: W = I.
%   tol   - (optional) tolerance vector for determining the order of 
%           reduced system, of the form [tol1, tol2], where:
%             tol1 specifies the tolerance for model reduction.
%                  Default: tol1 = NS*epsilon_machine*HSV(1), where NS
%                  is the order of the alpha-stable part of G.
%             tol2 specifies the tolerance for minimal realization.
%                  Default: tol2 = NS*epsilon_machine*HSV(1).
%   discr - (optional) type of system:
%             = 0 : continuous-time (default);
%             = 1 : discrete-time.
%   ord   - (optional) desired order of reduced system.
%             Default: ord = -1 (order determined automatically).
%   alpha - (optional) vector of the form [sdeg, alphac, alphao], where
%             sdeg is the stability boundary for the eigenvalues of A;
%               default:    -sqrt(epsilon_machine)  for continuous-time;
%                        1.0-sqrt(epsilon_machine)  for discrete-time;
%             alphac is the weighting factor for the frequency-weighted 
%                    controllability Grammian; default: alphac = 0;
%             alphao is the weighting factor for the frequency-weighted 
%                    observability Grammian;   default: alphao = 0.
%             Note: For further explanation of alphac and alphao, see
%                   subroutine AB09ID.
%
% Description of output parameters:
%   Ar, Br, 
%   Cr, Dr - matrices of the reduced system. 
%   HSV    - frequency-weighted Hankel singular values of the 
%            alpha-stable part.
%   info   - warning message code:
%            info = 1 - selected order greater than the order
%                       of a minimal realization;
%            info = 2 - selected order corresponds to repeated singular
%                       values, which are neither all included nor all 
%                       excluded from the reduced model;
%            info = 3 - selected order less than the order of
%                       the unstable part.
%            info = 10+K - K violations of the numerical stablity
%                       condition appeared during eigenvalue assignment.
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
