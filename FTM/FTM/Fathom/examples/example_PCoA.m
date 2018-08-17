% Examples for PCoA
% 
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: updated for new plotting routine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/legendre_gallagher.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is an artificial data set comprising 19 sample sites and 9 species:
% Y       = abundance values of species for each site
% yLabels = cell array of species names
% sLabels = cell array of site names
% 
% (Source: Figure 3 in LLegendre, P., and E. E. Gallagher. 2001.
% Ecologically meaningful transformations for ordination of species data. Oecologia 129: 271?280.) 

% Load the data:
load legendre_gallagher.mat;

% Create a Bray-Curtis dissimilarity matrix among observations:
dis = f_dis(Y,'bc');

% Perform Principal Coordinates Analysis, show diagnostic plots:
pcoa = f_pcoa(dis,1);

% Create ordination diagram and biplot vectors (use 'weighted' vectors
% because the data represents species abundance data):
f_pcoaPlot(pcoa,sLabels,Y,yLabels,0,'none',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/spiders.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% -----NOTES:-----
% This is data from van der Aart & Smeenk-Enserink (1975) and is Hunting spider
% abundances for 12 species (spiders) taken from 28 sites (labels) and
% associated environmental data (env).
% 
% [van der Aart, P. J. M. & N. Smeenk-Enserink. 1975. Correlations between
% distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
% characteristics in a dune area. Netherlands Journal of Zoology 25: 1?45.]  

%%%%%%%%%%%%%%%%%%%%%%
%    BIOTIC DATA:    %
%%%%%%%%%%%%%%%%%%%%%%
% Load the data files:
load spiders.mat;
 
% Create a Bray-Curtis dissimilarity matrix among observations:
dis = f_dis(spiders,'bc');
 
% Perform Principal Coordinates Analysis, show diagnostic plots:
pcoa = f_pcoa(dis,1);

% Create ordination diagram and biplot vectors (use 'weighted' vectors
% because the data represents species abundance data):
[h,vec] = f_pcoaPlot(pcoa,site_labels,spiders,spiders_labels,0,'none',1);

% Plot vectors for the 3 most important taxa for each axis:
[null,idx_I]  = sortrows(abs(vec),-1); % sort rows descending by column 1
[null,idx_II] = sortrows(abs(vec),-2); % sort rows descending by column 2
idx = f_unique([idx_I(1:3);idx_II(1:3)]);
idx = idx(1:3);
f_pcoaPlot(pcoa,site_labels,spiders(:,idx),spiders_labels(idx),0,'none',1);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ENVIRONMENTAL DATA:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a Euclidean distance matrix among observations:
eDis = f_dis(env,'euc');

% Perform Principal Coordinates Analysis, show diagnostic plots:
ePcoa = f_pcoa(eDis,1);

% Create ordination diagram and standard correlation vectors:
f_pcoaPlot(ePcoa,site_labels,env,env_labels,0,'none',0)
