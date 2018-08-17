% Examples for DISSIMILARITIES
% 
% by David L. Jones, Mar-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/spiders.mat'

% -----NOTES:-----
% This is data from van der Aart & Smeenk-Enserink (1975) and is Hunting spider
% abundances for 12 species (spiders) taken from 28 sites (labels) and
% associated environmental data (env).
% 
% [van der Aart, P. J. M. & N. Smeenk-Enserink. 1975. Correlations between
% distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
% characteristics in a dune area. Netherlands Journal of Zoology 25: 1?45.]  

% -----ANALYSIS:-----
% Load the data files:
load spiders.mat
 
% Create a Bray-Curtis dissimilarity matrix among observations:
disBC = f_dis(spiders,'bc');

% Jaccard index:
disJAC = f_dis(spiders,'jac');

% Compare the 2 distance measures:
r = f_mantel(disBC,disJAC,1)
% r = 0.9165
