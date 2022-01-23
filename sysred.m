%SYSRED  MEX-function based on SLICOT model reduction routines.
%
%   [Ar,Br,Cr,Dr,HSV] = SYSRED(METH,A,B,C,D,TOL,DISCR,ORDER,ALPHA) 
%
%   SYSRED returns for an original continuous- or discrete-time
%   state-space system (A,B,C,D) a reduced order state-space system
%   (Ar,Br,Cr,Dr) and the Hankel singular values HSV of the 
%   ALPHA-stable part. The order of the reduced model is determined 
%   either by the number of Hankel-singular values greater than TOL or
%   by the desired order ORDER. 
%   
%   Description of other input parameters: 
%   METH  - method flag with decimal form c*10+m, where:
%             m specifies the basic model reduction method;
%             c specifies the comprime factorization approach to be 
%               used in conjunction with the method specified by m.
%           Allowed values for m:
%             m = 1 : Balance & Truncate method with balancing
%             m = 2 : Balance & Truncate method (no balancing)
%             m = 3 : Singular Perturbation Approximation with balancing
%             m = 4 : Singular Perturbation Approximation (no balancing)
%             m = 5 : Optimal Hankel-Norm Approximation 
%           Allowed values for c (only for m = 1..4):
%             c = 0 : no coprime factorization is used (default)
%             c = 1 : RCF with inner denominator    
%             c = 2 : LCF with inner denominator    
%             c = 3 : RCF with ALPHA stability degree  
%             c = 4 : LCF with ALPHA stability degree  
%   TOL   - (optional) tolerance vector for determining the order of 
%           reduced system of form TOL = [tol1, tol2, tol3], where
%             tol1 specifies the tolerance for model reduction;
%                  default: tol1 = epsilon_machine*Hankel_norm(A,B,C)
%             tol2 specifies the tolerance for minimal realization in
%                  case of m = 3, 4 or 5;
%                  default: tol2 = epsilon_machine*Hankel_norm(A,B,C)
%             tol3 specifies the controllability/observability tolerance
%                  for computing coprime factorizations, as follows: 
%                  controllability tolerance in case c = 1 or 3;
%                  default: epsilon_machine*max(norm(A),norm(B));
%                  observability tolerance in case c = 2 or 4 
%                  default: epsilon_machine*max(norm(A),norm(C))
%   ORDER - (optional) desired order of reduced system
%             default: ORDER = -1 (order determined automatically)
%   DISCR - (optional) type of system
%              = 0 : continuous-time (default)
%              = 1 : discrete-time
%   ALPHA - (optional) stability boundary for the eigenvalues of A
%             default:    -sqrt(epsilon_machine)  for continuous-time
%                      1.0-sqrt(epsilon_machine)  for discrete-time

%   RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%	A. Varga 2-15-99.

%       Reference:
%       [1] A. Varga
%           Selection of Model Reduction Routines.
%           SLICOT Working Note 1998-2 
%           ftp://wgs.esat.kuleuven.ac.be/pub/WGS/REPORTS/SLWN1998-2.ps.Z
%
