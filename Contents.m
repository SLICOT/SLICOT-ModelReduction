% SLMODRED - SLICOT Model and Controller Reduction Toolbox.
% RELEASE 2.0 of SLICOT Toolboxes.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% This toolbox includes M- and MEX-functions for model and 
% controller reduction for standard state space systems.
%
% MATLAB M-files (except help files).
% -----------------------------------
%
% Model and controller reduction for standard state space systems.
% 
%   bta         - Balance & Truncate approximation without balancing. 
%   btabal      - Balance & Truncate approximation with balancing. 
%   bta_cf      - Balance and truncate approximation of coprime factors 
%                 without balancing. 
%   btabal_cf   - Balance and truncate approximation of coprime factors
%                 with balancing. 
%   hna         - Hankel-norm approximation.
%   spa         - Singular perturbation approximation without balancing.
%   spabal      - Singular perturbation approximation with balancing. 
%   spa_cf      - Singular perturbation approximation of coprime factors 
%                 without balancing.
%   spabal_cf   - Singular perturbation approximation of coprime factors 
%                 with balancing. 
%   bst         - Balanced stochastic truncation approximation. 
%   cfconred    - Coprime factorization based reduction of state feedback 
%                 controllers. 
%   fwbconred   - Frequency-weighted balancing related controller reduction. 
%   fwbred      - Frequency-weighted balancing related model reduction. 
%   fwhna       - Frequency-weighted Hankel-norm approximation. 
%   sysredget   - Get SYSRED OPTIONS parameters. 
%   sysredset   - Create/alter SYSRED OPTIONS structure. 
%
% MATLAB MEX-files.
% -----------------
%
% Model reduction for standard state space systems.
%   
%   sysred      - Computes for an original continuous- or discrete-time
%                 state-space system (A,B,C,D) a reduced order state-space system
%                 (Ar,Br,Cr,Dr) and the Hankel singular values of the 
%                 ALPHA-stable part.  Various techniques can be used.
%
% Model/Controller reduction for standard state space systems.
%
%   bstred      - Computes a reduced order state-space system using a 
%                 balanced stochastic truncation model reduction method. 
%   conred      - Computes a reduced order controller using a frequency-weighted 
%                 balancing related controller reduction method.
%   fwehna      - Computes a reduced order state-space system using a 
%                 frequency-weighted Hankel-norm approximation method.
%   fwered      - Computes a reduced order state-space system using a 
%                 frequency-weighted balancing related model reduction method.
%   sfored      - Computes a reduced order controller using a 
%                 state-feedback-observer controller reduction method.
%   
% Auxiliary MEX-file.
%
%   arebench    - Generates benchmark examples for algebraic Riccati equations
%                 (used in the demonstration file cmp_red). 
%
% Demonstrations.
% ---------------
%
%   slredemo    - Model and Controller Reduction Demo.
%   cmp_red     - Comparison between SLMODRED function hna and
%                 MATLAB Robust Control function 
%                 ohklmr   (for MATLAB version <  7) or
%                 hankelmr (for MATLAB version >= 7).
%
% Note:
% -----
%
%   Command-line help functions are available for all MATLAB M- and MEX-functions 
%   included in this toolbox.
%

%  ..CONTRIBUTOR..
%
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2004.
%
%   Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2005.
%   V. Sima, Jan. 2022.




