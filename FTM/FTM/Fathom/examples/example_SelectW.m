% Example of selecting among competing spatial weighting matrices:
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/oribatid_mites.mat'

% -----Notes:-----
% This example reproduces that provided for the 'test.W' function in 'Dray, S.
% 2006. Moran's eigenvectors of spatial weighting matrices in R', which is
% included in the 'spacemakeR for R' documentation. Available from:
% http://biomserv.univ-lyon1.fr/~dray/software.php

% Load data:
load oribatid_mites.mat

% Hellinger transform the species data:
bio.H = f_hellinger(bio.num);

% Perform preliminary RDA on X,Y coordinates to identify linear trends:
result_XY = f_rda(bio.H,[env.x env.y],[0],1000,1);

% Residuals serve as the de-trended data:
bio.dt = result_XY.res;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Test a binary (b) spatial weighting matrix:                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Delaunay Triangulation based neighbor graph:
tri = f_delaunay([env.x env.y]);

% Create a weighting matrix ('A') that doesn't change 'B':
A = ones(size(tri.B));

% Create Eigenvector Maps, no iterations, KEEP neg eigenvalues:
MEM_b = f_eigenMaps(tri,A,0,1);

% AIC-based stepwise selection of MEM's:
[MEM_b.model,MEM_b.cond] = f_rdaAIC(bio.dt,MEM_b.evects,0,0);

% Show the results:
MEM_b.model
% ans =
%        RSS: [69x1 double]
%         R2: [69x1 double]
%      R2adj: [69x1 double]
%        AIC: [69x1 double]
%        var: {1x69 cell}
%        idx: [3 2 4 7 1 5 16 53]
%     minAIC: -95.365
%       null: -87.471
%          X: [70x8 double]

% Results from 'spacemakeR:
% $best$AICc  = -90.922628 -92.465108 -93.783284 -94.725071 -95.195222
% $best$ord   =  3 2 7 1 16
% $best$AICc0 = -87.47112



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Test a spatial weighting matrix based on a linear function (f1):     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Delaunay Triangulation based neighbor graph:
tri = f_delaunay([env.x env.y]);

% Create a weighting matrix ('A') based on a linear distance function:
A = f_dis2sim(tri.dis,1);

% Create Eigenvector Maps, no iterations, KEEP neg eigenvalues:
MEM_f1 = f_eigenMaps(tri,A,0,1);

% AIC-based stepwise selection of MEM's:
[MEM_f1.model,MEM_f1.cond] = f_rdaAIC(bio.dt,MEM_f1.evects,0,0);

% Show the results:
MEM_f1.model
% ans =
%        RSS: [69x1 double]
%         R2: [69x1 double]
%      R2adj: [69x1 double]
%        AIC: [69x1 double]
%        var: {1x69 cell}
%        idx: [3 2 7 1 16 5 53]
%     minAIC: -95.801
%       null: -87.471
%          X: [70x7 double]
%
%

% Results from 'spacemakeR:
% $best$AICc  = -91.488003 -94.917598 -95.489469 -95.905881 -96.339367 -96.716195
% $best$ord   =  3 1 2 15 7 24
% $best$AICc0 = -87.47112



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Test a range of concave down weighting functions:                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Delaunay Triangulation based neighbor graph:
tri = f_delaunay([env.x env.y]);

% Initialize:
clear A;
e = [2:10]'; % range of exponents

% Create neighbor graphs and corresponding weighting matrices:
for i = 1:size(e,1)
   N{i} = tri;                       % neighbor graph
   A{i} = f_dis2sim(tri.dis,2,e(i)); % weighting matrix ('A')
end

% Select optimal spatial weighting matrix, keep negative eigenvalues:
model_f2 = f_selectW(bio.dt,N,A,1);
% model_f2 = 
%        RSS: [69x1 double]
%         R2: [69x1 double]
%      R2adj: [69x1 double]
%        AIC: [69x1 double]
%        var: {1x69 cell}
%        idx: [3 2 7 1 16 4 5 53]
%     minAIC: -95.653
%       null: -87.471
%          X: [70x8 double]
clear i e;

% Results from 'spacemakeR:
% Best model:
% 
% 
%    y dmax      AICc NbVar
% 9 10    3 -96.22122     6
% 
% $all
%    y dmax      AICc NbVar
% 1  2    3 -95.90200     5
% 2  3    3 -95.68185     7
% 3  4    3 -95.69664     6
% 4  5    3 -95.88353     6
% 5  6    3 -96.02857     6
% 6  7    3 -96.12361     6
% 7  8    3 -96.17905     6
% 8  9    3 -96.20812     6
% 9 10    3 -96.22122     6


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Final results:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% For the 9 concave-down spatial weighting functions (f2), e=2 was determined to
% provide the optimal spatial model; although, 'spacemakeR' ranked e=10 as the
% best. However, BOTH platforms ranked 'f1' as the best among all 11 competing
% models (i.e., b, f1, & the f2's).
