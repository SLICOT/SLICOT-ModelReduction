echo on,
%   This demonstration example uses the data for the 16-th order demo
%   example from MATLAB Robust Control Toolbox.  This system has
%   two inputs and two outputs.  The A matrix is block-diagonal 
%   with two 8-by-8 diagonal blocks.
 
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

%       Print the model:

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
a
echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
b
echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
c
echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
d
echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Find a 6-th reduced order model using SLICOT 
%       optimal Hankel norm approximation interface to MATLAB.
%

[ sysr, slhsv ] = hna( sys, [], 6 );

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Print and plot the Hankel singular values:

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp( 'Hankel singular values' )
disp( slhsv )
echo on
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

disp( ' ' )
disp( [ 'H-infty error norm for the 6-th order approximation = ', ...
        num2str( 2*sum( slhsv(7:16) ) ) ] )
echo on

%       Singular values plot.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off

[ svo, w ] = sigma( sys ); svo = 20*log10( svo );
svr = sigma( sysr, w );    svr = 20*log10( svr );
h = axes;
semilogx( w, svo(1,:), 'r', w, svo(2,:), 'r--', ...
          w, svr(1,:), 'b', w, svr(2,:), 'b--' )
set( h, 'FontSize', 12 )
title( 'Optimal Hankel Norm Approximation (order 16 --> order 6)' )
xlabel( 'Frequency (rad/sec)' )
ylabel( 'Singular values (dB)' )
legend( h, 'Original system ({\it n} = 16), y_1', 'idem, y_2', ...
           'Reduced system ({\it n_r} = 6), y_1', 'idem, y_2', 'location', 'Best' )
%grid on
shg

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off
close( gcf )
echo on

%       Similarly, find a 4-th reduced order model.
%

[ sysr4, slhsv ] = hna( sys, [], 4 );

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

echo off
disp( ' ' )
disp( [ 'H-infty error norm for the 4-th order approximation = ', ...
         num2str( 2*sum( slhsv(5:16) ) ) ] )
echo on

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Singular values plot.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off

svr4 = sigma( sysr4, w );  svr4 = 20*log10( svr4 );
h = axes;
semilogx( w, svo(1,:),  'r', w, svo(2,:),  'r--', ...
          w, svr4(1,:), 'b', w, svr4(2,:), 'b--' )
set( h, 'FontSize', 12 )
title( 'Optimal Hankel Norm Approximation (order 16 --> order 4)' )
xlabel( 'Frequency (rad/sec)' )
ylabel( 'Singular values (dB)' )
legend( h, 'Original system ({\it n} = 16), y_1', 'idem, y_2', ...
           'Reduced system ({\it n_r} = 4), y_1', 'idem, y_2', 'location', 'Best' )
shg

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off
close( gcf )
echo on

%       Use balanced stochastic truncation with BstBeta = 1.
%

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

oldoptions = sysredset( 'Order', 6, 'bstb', 1 );
[ sysb, hsvb ] = bst( sys, oldoptions );

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Additionally, use singular perturbation approximation
%       (BalredMethod = SPA).
%

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

options = sysredset( oldoptions, 'balr', 'spa' );
[ sysbs, hsvbs ] = bst( sys, options );

if pause_wait < 0,  pause,  else  pause(pause_wait),  end

%       Singular values plot.

if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off

svb = sigma( sysb,  w );  svb = 20*log10( svb );
svs = sigma( sysbs, w );  svs = 20*log10( svs );
h = axes;
semilogx( w, svo(1,:), 'r', w, svo(2,:), 'r--', ...
          w, svr(1,:), 'b', w, svr(2,:), 'b--', ...
          w, svb(1,:), 'g', w, svb(2,:), 'g--', ...
          w, svs(1,:), 'k', w, svs(2,:), 'k--' )
set( h, 'FontSize', 12 )
title( 'Several order reduction methods (order 16 --> order 6)' )
xlabel( 'Frequency (rad/sec)' )
ylabel( 'Singular values (dB)' )
legend( h, 'Original system, y_1', 'idem, y_2', ...
           'HNA, y_1', 'idem, y_2', ...
           'BST (BstBeta = 1), y_1', 'idem, y_2', ...
           'BST (+ Balr = SPA), y_1', 'idem, y_2', 'location', 'Best' )
shg

echo on
if pause_wait < 0,  pause,  else  pause(pause_wait),  end
echo off
close( gcf )
