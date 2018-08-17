% Example of using Moran's eigenvector maps (MEM's) in spatial ecology
% 
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/oribatid_mites.mat''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Oribatid Mite Example of Griffith & Peres-Neto, 2006:            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Load data:
load oribatid_mites.mat

% Hellinger transform the species data:
bio.H = f_hellinger(bio.num);

% Dummy code qualitative variables:
env.Qsubstra = f_dummy(env.substra,1);
env.Qschrub  = f_dummy(env.schrub,1);
env.Qmicro   = f_dummy(env.micro,1);

% Perform preliminary RDA on X,Y coordinates to identify linear trends:
result_XY = f_rda(bio.H,[env.x env.y],[0],1000,1);

% Residuals serve as the de-trended data:
bio.dt = result_XY.res;

% Create a Minimum Spanning Tree:
mst = f_mst(f_euclid([env.x env.y]'),[env.x env.y]);

% Create a weighting matrix ('A') based on the truncated distance matrix:
A = f_dis2sim(mst.tDis,2,2);

% Create Eigenvector Maps, use 1000 iterations, don't keep neg eigenvalues:
MEM = f_eigenMaps(mst,A,1000,0);

% Set up environmental variables:
X = [env.density env.water env.Qsubstra env.Qschrub env.Qmicro];

% Use Pedro's routine to select MEM's:
% pedro = StepWiseMoran(bio.dt,MEM.evects,MEM.W,X,1000,0.05);


% Select MEM's that contribute to explaining the spatial structure:
f_eigenMapsStepwise(bio.dt,X,MEM,1000,0.99,1);