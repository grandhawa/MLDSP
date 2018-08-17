% Compare results obtained by Matlab vs. R
% 
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
 
% File: '.../examples/Matlab_v_R.mat'

% Load the data file:
load Matlab_v_R;

% Compare response data:
sum(sum(bio.dt - faudt))
% ans = 1.1429e-14
% 
f_procrustes(bio.dt,faudt,0,0,0)
% ans = 4.4409e-16

% Compare MEM's:
sum(sum(MEM_f1.evects - MEM))
% ans = -7.051e-14
% 
f_procrustes(MEM_f1.evects,MEM,0,0,0)
% ans = 2.2204e-15

% Note: differences between the two platforms are VERY slight, which in turn
% results in slightly different results when comparing the results of model
% selection in f_rdaAIC vs. ortho.aic/test.W
