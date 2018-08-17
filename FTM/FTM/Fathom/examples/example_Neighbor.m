% Example of creating Neighborhood Graphs
% 
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/neighbors.mat'

% -----Notes:-----
% This example is similar to that presented in 'Dray, S. 2006. Moran's
% eigenvectors of spatial weighting matrices in R', which is included in the
% 'spacemakeR for R' documentation. Available from:
% http://biomserv.univ-lyon1.fr/~dray/software.php  

% Load the data file:
load neighbors.mat

% Create Neighborhood Graphs:
tri  = f_delaunay(dat);
gab  = f_gabriel(dat);
rel  = f_relNeigh(dat);
mst  = f_mst(f_euclid(dat'),dat);
dnn1 = f_dnn(f_euclid(dat'),[0 0.3],dat);
dnn2 = f_dnn(f_euclid(dat'),[0 0.5],dat);

% Plot the graphs:
figure;

subplot(3,3,1)
f_plotNeigh(tri,1,'r');
title('\bfDelaunay Triangulation');

subplot(3,3,2)
f_plotNeigh(gab,1,'r');
title('\bfGabriel Graph');

subplot(3,3,3)
f_plotNeigh(rel,1,'r');
title('\bfRelative Neighbor Graph');

subplot(3,3,4)
f_plotNeigh(mst,1,'r');
title('\bfMinimum Spanning Tree');

subplot(3,3,5)
f_plotNeigh(dnn1,1,'r');
title('\bfNeighbors if 0<d<0.3');

subplot(3,3,6)
f_plotNeigh(dnn2,1,'r');
title('\bfNeighbors if 0<d<0.5');
