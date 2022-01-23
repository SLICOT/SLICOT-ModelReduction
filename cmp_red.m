% Script for comparing model reduction functions for continuous- and
% discrete-time models from the AREs benchmark collection, using SLICOT 
% and MATLAB solvers implemented in the function red_slv.
%
% Existing variables are saved, then cleared, and finally restored.

%   RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2004.
%
%   Revisions:
%   V. Sima, Apr. 2005.
%   P. Vacher, July 2005.
%   V. Sima, Nov. 2005, Mar. 2009.
%

% Save and then clear variables.

vars  = who( '*' );
saved = 0;

if ~isempty( vars ),
    clear vars
    if exist( 'currwork.mat', 'file' ) ~= 2,
        save currwork
        saved = 1;  save savevar saved
    else
        disp( ' ' )
        disp( [ 'Current directory is ', pwd ] )
        disp( ' ' )
        disp( 'Variables will be saved, then cleared, and finally restored.' );
        matfile = input( 'Specify a MAT-file name for saving the workspace, MAT_name = ', 's' );
        save( matfile )
        saved = 2;  save savevar saved matfile
    end
else
    save savevar saved
end

clear variables;

% Changing directory.

curdir = pwd;
cd( fileparts( which( 'cmp_red' ) ) );

% Set the example type flag Ts:
% Ts = 0, for continuous-time systems.
% Ts = 1, for discrete-time systems.

% Continuous-time systems

Ts = 0;

space = [];  l = 0;

% Set the option for plotting the Hankel singular values,
% in order to choose the reduced system order, nr.  
% Default:  nr is already set in this script for each system.

if ~exist( 'setnr', 'var' ) || isempty( setnr ),  setnr = 0;  end

% Set solver activation flag.  All are used, by default.  
% To select some solver(s), set UseAll = 0, and the corresponding
% flag(s) to 1.

UseAll    = 1;
UseSLhna  = 0;
UseMThna  = 0;
% UseAll    = 0;

% Define the MATLAB function to be used, depending on MATLAB version.

vers7 = version;  vers7 = strcmp( '7', vers7(1) );

% Fictitious calls for loading all solvers. This avoids counting the
% loading time for timing the first run.

disp( 'Fictitious calls for loading the solvers.' )
disp(' ')

A = rand( 2,2 );   B  = rand(2,1);  C = rand(1,2);  D = rand(1,1) ;
n = size( A, 1 );  nr = n;

sys = ss( A, B, C, D );

[ sysr, slhsv ]   = hna( sys, [], nr );
if vers7,
    [ sysm, redinfo ] = hankelmr( sys, nr );        nmfct = 'hankelmr';
else
    [ sysm, boundM, mhsv ] = ohklmr( sys, 1, nr );  nmfct = 'ohklmr';
end

disp(' ')
disp('Comparison between SLICOT and MATLAB model reduction solvers')
disp('------------------------------------------------------------')
disp(' ')

[X,A,R,W,B,C] = arebench(1,1,6,[0 0]);  nr1 = 1;  nr2 = 6;  nr = 24;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;  

[X,A,R,W,B,C] = arebench(1,2,9,[0 0]);  nr1 = 2;  nr2 = 9;  nr = 49;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

[X,A,R,W,B,C] = arebench(1,3,1,[0 0],[],5);  nr1 = 3;  nr2 = 1;  nr = 8;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

[X,A,R,W,B,C] = arebench(1,4,2,[0 0]);  nr1 = 4;  nr2 = 2;  nr = 22;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

[X,A,R,W,B,C] = arebench(1,4,2,[0 0],[.05 .1 .1 .1 .5 .1 .5],20);  nr = 14;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

[X,A,R,W,B,C] = arebench(1,4,4,[0 0]);  nr1 = 4;  nr2 = 4;  nr = 300;  
disp( [ sprintf( '\n' ), 'CARE example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

% Discrete-time systems

Ts = 1;

[X,A,R,W,B,C] = arebench(2,1,11,[0 0 0]);  nr1 = 1;  nr2 = 11;  nr = 8; 
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

[X,A,R,W,B,C] = arebench(2,1,12,[0 0 0]);  nr1 = 1;  nr2 = 12;  nr = 8;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

[X,A,R,W,B,C] = arebench(2,1,13,[0 0 0]);  nr1 = 1;  nr2 = 13;  nr = 16;
disp( [ sprintf( '\n' ), 'Example ', num2str( nr1 ), '.', num2str( nr2 ),'.' ] )
red_slv;

disp( ' ' )
disp( '--------------------------------------------------')
disp( '                                       Time       ')   
disp( '--------------------------------------------------')
disp( [ '  nr1 nr2   n    m    p   nr   Ts   hna  ', nmfct ] )
disp( '--------------------------------------------------')
		                    
disp( [space, int2str( prob(:,1) ), space, int2str( prob(:,2) ), ...
       space, int2str( dims(:,1) ), space, int2str( dims(:,2) ), ... 
       space, int2str( dims(:,3) ), space, int2str( dims(:,4) ), ...
       space, int2str( dims(:,5) ), ...
       space, num2str( timing(:,1), '%5.2f' ), ...
       space, num2str( timing(:,2), '%5.2f' ) ] )	
disp( '--------------------------------------------------')
		                    
disp( ' ' )
disp( '-------------------------------------------------------------------------------')
disp( '                                                  Error bounds                 ')   
disp( '-------------------------------------------------------------------------------')
disp( [ '  nr1 nr2   n    m    p   nr   Ts    SLICOT hna  MATLAB ', nmfct, ' RelErr_S2M ' ] )
disp( '-------------------------------------------------------------------------------')
		                    
disp( [space, int2str( prob(:,1) ),   space, int2str( prob(:,2) ), ...
       space, int2str( dims(:,1) ),   space, int2str( dims(:,2) ), ... 
       space, int2str( dims(:,3) ),   space, int2str( dims(:,4) ), ...
       space, int2str( dims(:,5) ),   ...
       space, num2str( err(:,1), 5 ), space, num2str( err(:,2), 5 ), ... 
       space, num2str( err(:,3), 5 ) ] )	
disp( '-------------------------------------------------------------------------------')
disp( ' ')		                    		                    

disp(   'Note: RelErr_S2M denotes the relative error (in the Euclidean norm)' )
disp(   '      between the Hankel singular values computed by SLICOT hna' )
disp( [ '      and MATLAB ', nmfct,'.' ] )

% Plot speed-ups of SLICOT over MATLAB.

disp( ' ')		                    		                    
disp( 'Plot speed-ups of SLICOT over MATLAB' )

tmp  = timing(:,1);
indx = find( tmp == 0 );
tmp(indx) = 0.01;            % Zero denominators are set to 0.01 sec.
speedup = timing(:,2)./tmp;

close( gcf )
set( axes, 'FontSize', 16 )
bar( speedup )
title( [ 'Speed-up factors of SLICOT hna over ', nmfct ] )
shg

disp( 'Press any key to continue' )
disp( ' ' )
pause
close( gcf )

% Restoring initial directory.

cd( curdir )

clear variables

% Load the initial variables, if needed. 

if exist( 'savevar.mat', 'file' ) == 2,
    load savevar;  delete( 'savevar.mat' );
    if saved == 1,
        load currwork;   delete( 'currwork.mat' );
    elseif saved == 2,
        load( matfile );  delete( [ matfile, '.mat' ] );  clear matfile;
    end
    clear saved
end
