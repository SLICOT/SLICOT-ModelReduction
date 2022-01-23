%   The SLICOT Model and Controller Reduction Toolbox contains 
%   tools for finding reduced-order models of linear multivariable 
%   dynamical systems.
%
%   SLICOT Model and Controller Reduction Toolbox demonstrations:
%   1) A detailed example: a 16-th order demo example.
%   2) Comparison of SLICOT and MATLAB optimal Hankel norm approximation
%      model reduction functions for the 16-th order demo example.
%   3) Comparison between SLICOT Model and Controller Reduction
%      Toolbox function hna and MATLAB Robust Control functions:
%        * ohklmr   for MATLAB versions <  7;
%        * hankelmr for MATLAB versions >= 7.
%
%   0) Quit.
%
%   Note: SLICOT Model and Controller Reduction Toolbox uses the system object 
%         from Matlab Control System Toolbox.
%

%   RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   V. Sima, Research Institute for Informatics, Bucharest, Nov. 2004.
%
%   Revisions:
%   P. Vacher 01-07-2005.
%   V. Sima 03-07-2020.
%

global pause_wait  % This could be used in pause(n) command.
                   % If pause_wait < 0, standard command pause is used (default).
                   % Any key should then be pressed to continue.
                   % It may need to use the command "pause on" before
                   % calling fstdemo.

k = 0;
while (1)
   disp(' ')
   help slredemo
   k = input('Select a demonstration example number: ');
   disp(' ')
   if isempty(k), k = 20; end
   if k == 0,  break,     end
   if k == 1,  slredemo1,  end
   if k == 2,  slredemo2,  end
   if k == 3,  cmp_red,    end
end
