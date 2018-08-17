% Example of Biotic - Environmental Correlation (aka BIOENV or BEST)
% 
% by David L. Jones, Jan-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% updated: Jun-2012
 
% File: '.../examples/spiders.mat'

% These data are from van der Aart & Smeenk-Enserink (1975) and are Hunting
% spider abundances for 12 species (spiders) taken from 28 sites
% (site_labels) and associated environmental data (env).
% 
% van der Aart, P. J. M. & N. Smeenk-Enserink. 1975. Correlations between
% distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
% characteristics in a dune area. Netherlands Journal of Zoology 25: 1-45.

load spiders.mat;                         % load the data
spiders2        = f_transform(spiders,3); % log-transform abundances
dis             = f_dis(spiders2,'bc');   % Bray-Curtis dissimilarity matrix

[res,resLabels] = f_bioenv(dis,env,env_labels,0,0,5,1);

% There are 63 possible subsets of 6 variables 
%   Processing    6 subsets of   1 variables 
%   Processing   15 subsets of   2 variables 
%   Processing   20 subsets of   3 variables 
%   Processing   15 subsets of   4 variables 
%   Processing    6 subsets of   5 variables 
%   Processing    1 subsets of   6 variables 
% 
%  ========================================== 
%  Rho    Variables 
%  ========================================== 
% 
%  1
%  0.7029 water
%  0.5454 light
%  0.5217 sand
%  0.4977 cover
%  0.4210 twigs
%  0.2493 herbs
% 
%  2
%  0.7527 twigs water
%  0.7484 light water
%  0.7096 light sand
%  0.7006 cover water
%  0.6977 sand twigs
% 
%  3
%  0.8030 sand twigs water
%  0.7755 light sand water
%  0.7572 cover sand twigs
%  0.7566 cover twigs water
%  0.7468 cover light sand
% 
%  4
%  0.8186 cover sand twigs water
%  0.7974 light sand twigs water
%  0.7884 cover light sand water
%  0.7829 cover herbs sand water
%  0.7727 herbs light sand water
% 
%  5
%  0.8092 cover light sand twigs water
%  0.8067 cover herbs sand twigs water
%  0.8017 cover herbs light sand water
%  0.7758 herbs light sand twigs water
%  0.7375 cover herbs light twigs water
% 
%  6
%  0.8010 cover herbs light sand twigs water
