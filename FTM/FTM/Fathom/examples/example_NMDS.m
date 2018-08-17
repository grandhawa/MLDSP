% Examples for nonmetric Multidimensional Scaling
% 
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% File: '.../examples/spiders.mat'

% -----NOTES:-----
% This is data from van der Aart & Smeenk-Enserink (1975) and is Hunting spider
% abundances for 12 species (spiders) taken from 28 sites (labels) and
% associated environmental data (env).
% 
% [van der Aart, P. J. M. & N. Smeenk-Enserink. 1975. Correlations between
% distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
% characteristics in a dune area. Netherlands Journal of Zoology 25: 1-45.]  

% -----ANALYSIS:-----
% Load the data files:
load spiders.mat
 
% Create a Bray-Curtis dissimilarity matrix among observations:
dis = f_dis(spiders,'bc');
 
% 2-D nonmetric Multidimensional Scaling:
mds_2 = f_nmds(dis,2,1,1);
% > 13 iterations, Final stress criterion = 0.0825375
% 
% Create ordination diagram and biplot vectors (use 'weighted' vectors
% because the data represents species abundance data):
f_nmdsPlot(mds_2,site_labels,spiders,spiders_labels,0,'none',1);

% 1-D nonmetric Multidimensional Scaling, plot configuration:
mds_1 = f_nmds(dis,1,1,1);
% > 6 iterations, Final stress criterion = 0.192875
% 
% Create ordination diagram and biplot vectors (use 'weighted' vectors
% because the data represents species abundance data):
f_nmdsPlot(mds_1,site_labels,spiders,spiders_labels,0,'none',1);
