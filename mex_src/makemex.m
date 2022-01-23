function makemex
%
% Function for generating gateway functions from MEX files for SLICOT Model and Controller Reduction Toolbox
% (Linux version, 64 bit).
% This file should be located in the directory containing the source MEX files, which is a subdirectory
% of the directory where the libraries slicot.a is located.
%
flags = 'FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer -fdefault-integer-8" -largeArrayDims';
slmodred_mex_src = '';
libslicot = '../slicot.a';
slmodred_mex = {
    'arebench.f', ...
    'bstred.f', ... 
    'conred.f', ...
    'fwehna.f', ...
    'fwered.f', ...
    'sfored.f', ... 
    'sysred.f', ...
    };
%
for k = 1:length(slmodred_mex)
    file = slmodred_mex{k};
    fprintf( 'mex %s %s%s.f %s %s -lmwlapack -lmwblas\n', flags, slmodred_mex_src, file, libslicot );
    %eval( sprintf( 'mex %s %s%s.f %s %s -lmwlapack -lmwblas\n', flags, slmodred_mex_src, file, libslicot ) );
end
