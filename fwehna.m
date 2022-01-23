% FWEHNA.F - MEX-function for SLICOT model reduction routine
%            AB09JD.F. 
%
% Matlab call:  
%   [Ar,Br,Cr,Dr,HSV,info] = FWEHNA(meth,A,B,C,D,V,W,tol,discr,ord,alpha)
%
% Purpose:
%   To compute a reduced order model (Ar,Br,Cr,Dr) for an original
%   state-space representation (A,B,C,D) by using the frequency 
%   weighted optimal Hankel-norm approximation method. 
%   The Hankel norm of the weighted error 
%     
%         op(V)*(G-Gr)*op(W)   
%     
%   is minimized, where G and Gr are the transfer-function matrices
%   of the original and reduced systems, respectively, V and W
%   are invertible transfer-function matrices of the left and right 
%   frequency weights, and op(X) denotes X, conj(X), inv(X), or 
%   conj(inv(X)). V and W are specified by their state space 
%   realizations (AV,BV,CV,DV) and (AW,BW,CW,DW), respectively. 
%   When minimizing the weighted error V*(G-Gr)*W, V and W must have 
%   poles distinct from those of G.
%   When minimizing conj(V)*(G-Gr)*conj(W), conj(V) and conj(W) must 
%   have poles distinct from those of G. 
%   Additionally, V and W must be invertible transfer-function matrices.
%   If the original system is unstable, then the frequency-weighted
%   Hankel-norm approximation is computed only for the 
%   alpha-stable part of the system. 
%   The order of the reduced model is determined either by the number
%   of frequency-weighted Hankel-singular values HSV greater than tol or 
%   by the desired order ord.
%
% Description of input parameters: 
%   meth  - method flag of decimal form ijk, where:
%             i = 1 : op(V) = V;
%             i = 2 : op(V) = conj(V);
%             i = 3 : op(V) = inv(V);
%             i = 4 : op(V) = conj(inv(V));
%             j = 1 : op(W) = W;
%             j = 2 : op(W) = conj(W);
%             j = 3 : op(W) = inv(W);
%             j = 4 : op(W) = conj(inv(W));
%             k = 1 : if possible, use standard inverse based method;
%             k = 2 : use inverse free method;
%             k = 3 : use an automatic method.
%   A,B,
%   C,D   - state-space system matrices of size N-by-N, N-by-M, P-by-N,
%           and P-by-M, respectively.
%   V     - (optional) contains the system matrix [AV BV; CV DV] of 
%           the left weighting.
%   W     - (optional) contains the system matrix [AW BW; CW DW] of 
%           the right weighting.
%   tol   - (optional) tolerance vector for determining the order of 
%           reduced system, of the form [tol1, tol2], where:
%             tol1 specifies the tolerance for model reduction.
%                  Default: tol1 = NS*epsilon_machine*Hankel_norm(G1s),
%                  where Hankel_norm(G1s) is the Hankel-norm of the 
%                  projection G1s of op(V)*G1*op(W) (see AB09JD, 
%                  Section METHOD) and NS is the order of G1s. 
%             tol2 specifies the tolerance for computing a minimal
%                  realization of the alpha-stable part of the 
%                  weighted original system.
%                  Default: tol2 = NS*epsilon_machine*Hankel_norm(G1s).
%   discr - (optional) type of system:
%              = 0 : continuous-time (default);
%              = 1 : discrete-time.
%   ord   - (optional) desired order of reduced system.
%             Default: ord = -1 (order determined automatically).
%   alpha - (optional) stability boundary for the eigenvalues of A.
%             Default:    -sqrt(epsilon_machine)  for continuous-time;
%                      1.0-sqrt(epsilon_machine)  for discrete-time.
%
% Description of output parameters:
%   Ar, Br, 
%   Cr, Dr - state-space matrices of the reduced system Gr.
%   HSV    - Hankel singular values of the projection G1s of
%            op(V)*G1*op(W) (see AB09JD, Section METHOD), 
%            where G1 is the alpha-stable part of G. 
%   info   - warning message code:
%            info = 1 - selected order greater than the order
%                       of a minimal realization;
%            info = 2 - selected order less than the order of
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
