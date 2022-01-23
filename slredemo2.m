echo on,
%   This demonstration compares some SLICOT and MATLAB model reduction
%   functions using the data for the 16-th order demo example from
%   MATLAB Robust Control Toolbox.  This system has two inputs and
%   two outputs.  The A matrix is block-diagonal with two 8-by-8 diagonal 
%   blocks.
 
echo off
%   RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2004,
%            Feb. 2009.
%
%   Revisions:
%   V. Sima, Jan. 2018.
%
echo on
global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.

if ~exist('pause_wait', 'var') || isempty(pause_wait),  pause_wait = -1;  end

%       Load the input-output data:

load MRdemo2M

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

sys = ss( a, b, c, d );

%       Find a 6-th reduced order model using SLICOT 
%       optimal Hankel norm approximation interface to MATLAB.
%

time = cputime;
[ sysr, slhsv ] = hna( sys, [], 6 );
timeS = cputime - time;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Plot the Hankel singular values:

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
close( gcf )
set( axes, 'FontSize', 16 )
bar( slhsv )
title( 'Hankel singular values' )
shg

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off
close( gcf )
echo on

%       Find a 6-th reduced order model using MATLAB 
%       optimal Hankel norm approximation.
%

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

time = cputime;
[ sysm, bound, mhsv ] = ohklmr( sys, 1, 6 );
timeM = cputime - time;

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Singular values plot.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off

[ svo, w ] = sigma( sys );  svo = 20*log10( svo );
svr = sigma( sysr, w );     svr = 20*log10( svr );
svm = sigma( sysm, w );     svm = 20*log10( svm );
h = axes;
semilogx( w, svo(1,:), 'r', w, svo(2,:), 'r--', ...
          w, svr(1,:), 'b', w, svr(2,:), 'b--', ...
          w, svm(1,:), 'g', w, svm(2,:), 'g--' )
set( h, 'FontSize', 12 )
title( 'SLICOT and MATLAB optimal HNA (order 16 --> order 6)' )
xlabel( 'Frequency (rad/sec)' )
ylabel( 'Singular values (dB)' )
legend( h, 'Original system, y_1', 'idem, y_2', ...
           'SLICOT HNA, y_1', 'idem, y_2', ...
           'MATLAB HNA, y_1', 'idem, y_2', 'location', 'Best' )
shg

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off
close( gcf )

disp( ' ' )
disp( [ 'SLICOT H-infty error norm for the 6-th order approximation = ', ...
        num2str( 2*sum( slhsv(7:16) ) ) ] )
disp( [ 'MATLAB H-infty error norm for the 6-th order approximation = ', ...
        num2str( bound ) ] )
disp( ' ' )
disp( [ 'SLICOT execution time = ', num2str( timeS ) ] )
disp( [ 'MATLAB execution time = ', num2str( timeM ) ] )
