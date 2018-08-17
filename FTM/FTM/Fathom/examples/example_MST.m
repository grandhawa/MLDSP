% Examples for Minimum Spanning Trees
% 
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/spiders.mat'
 
% -----Notes:-----
% This is data is from from van der Aart & Smeenk-Enserink (1975) and is Hunting
% spider abundances for 12 species (spiders) taken from 28 sites (site_labels)
% and associated environmental data (env). 
 
% Load data file: 
load spiders.mat
% 
% Create dissimilarity matrix:
dis = f_dis(spiders,'bc');
% 
% Obtain a Euclidean-embedding via PCoA:
pcoa = f_pcoa(dis,0,0,0);

% Create Minimum Spanning Tree:
mst = f_mst(dis,pcoa.scores);

% Plot the tree:
f_plotNeigh(mst);
title('\bfMinimum Spanning Tree');
xlabel('Axis I');
ylabel('Axis II');
axis([-0.22 0.40 -0.30 0.33]);
