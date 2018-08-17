% Examples for Principal Component Analysis
% by David L. Jones, 2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% Jan-2013: updated

% File: '.../examples/spiders.mat'

% -----References:-----
% van der Aart, P. J. M. and N. Smeenk-Enserink. 1975. Correlations between
%   distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
%   characteristics in a dune area. Netherlands Journal of Zoology 25: 1-45.  

% Load the data file:
load spiders.mat

% PCA using correlation matrix:
pca = f_pca(env,1,2);

% Create plot:
f_pcaPlot(pca,[],env_labels,0.05);