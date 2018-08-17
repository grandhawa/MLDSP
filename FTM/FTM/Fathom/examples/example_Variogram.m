% Example of creating a multivariate empirical variogram
% 
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Note: this example follows the oribatid mite illustration (Fig. 3) of
% Dray et al., 2006

% Load the data file:
load variogram.mat;

% Create Bray-Curtis dissimilarity matrix:
yDis = f_dis(Y,'bc');

% Calculate variogram using Sturge's Rule, plot results:
result = f_variogram(yDis,X,0,1000,1);

% Show number of connections within each distance class:
disp('# Connections:'); fprintf('%d ',result.nc); disp(' ');
% # Connections:
% 66 70 70 70 70 70 70 70 65 48 40 31 17 9

% Show number of members within each distance class:
disp('# Members:'); fprintf('%d ',result.nw); disp(' ');
% # Members:
% 110 269 383 338 287 228 178 181 148 112 78 60 29 14  
