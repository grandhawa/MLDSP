% Examples for RDA
% 
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/legendre_113.mat'

% -----NOTES:-----
% This file reproduces the RDA analysis of coral reef fish abundances in Table
% 11.3 of Legendre & Legendre (1998); of particular interest is the use of the
% qualitative explanatory variable 'substrate type'.

% -----ANALYSIS:-----
% Load data file:
load legendre_113.mat

% Create categorical variable 'substrate type' (1=coral, 2=sand, 3=other)
sub    = [2 2 2 3 1 3 1 3 1 3]';
x_txt  = {'depth' 'coral' 'sand' 'other'}';
 
% Dummy code the substrate type (trim last column to avoid a singular matrix):
sub_dum = f_dummy(sub,1);

% Set up x:
x = [depth sub_dum];
 
% Perform RDA on abundances, don't use last column of dummy codes:  
result = f_rda(y,x,0,1000,1);
% 
% Permuting residuals under a full model 999 times...
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 47.6417    p    =  0.00100 
% R2 = 0.9597   R2adj =  0.93957 
% No. of permutations = 1000 
% --------------------------------------------------
% Canonical Eigenvalues:
%   74.5227  24.9420  8.8761
% Residual Eigenvalues:
%   4.1888  0.3139  0.0370  0.0085
% Species-Environment Correlations (r):
%   0.9986  0.9965  0.9797
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.9597): 
%   0.6601  0.2209  0.0786
% Cumulative: 
%   0.6601  0.8811  0.9597
% Residual axes  (total = 0.0403):
%   0.0371  0.0028  0.0003  0.0001
% Cumulative: 
%   0.0371  0.0399  0.0402  0.0403
% ==================================================


%-----RDA = biplot vectors for X,Y:-----
% Replace trimmed X with untrimmed version, standardize:
result.X = f_stnd([depth f_dummy(sub,0)]); 
% 
f_rdaPlot(result,y,0,10,0.25,y_txt,x_txt);


%-----RDA = biplot vectors for Depth,Y & centroids for Substrate:-----
% Replace X with just depth, standardize:
result.X = f_stnd(depth); 
% 
f_rdaPlot(result,y,[0],10,0.25,y_txt,x_txt(1));
% 
% Plot binary variables in X as centroids:
C = f_centroid(result.siteScores,sub);
h = text(C(:,1),C(:,2),x_txt(2:end));
   set(h,'FontSize',8,'HorizontalAlignment','center','Color','r');
clear C h;



%-----RDA = biplot vectors for X, WAScores for Y:-----
% Replace X with untrimmed version, standardize:
result.X = f_stnd([depth f_dummy(sub,0)]); 
% 
f_rdaPlot(result,y,[1],10,0.25,y_txt,x_txt);

