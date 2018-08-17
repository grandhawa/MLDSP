% Examples for Procrustes Analysis
% 
% by David L. Jones, 2003
% 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/losos.mat'

% This example follows example 1 of Peres-Neto & Jackson (2001) and uses
% data from Table 1 of Losos (1990).
% 
% Load the data file:
load losos.mat
% 
% Log transform the variables:
morph_log   = f_transform(morph,3);
perform_log = f_transform(perform,3);
% 
% PCA on correlation matrix:
morph   = f_pca(morph_log,0,2);
perform = f_pca(perform_log,0,2);
% 
% Scale variance of scores on each axis = 1:
morph.scores = f_stnd(morph.scores);
% 
% Procrustes Analysis:
[m2,morph_scl,perform_scl,resid,prob,H] = f_procrustes(morph.scores(:,1:2),...
   perform.scores(:,1:2),1,1000,1);
% 
% Permuting the data 999 times...
% 
prob
% 
% prob =   0.0010
% 
% Add labels to Superimposition plot:
f_labelplot([morph_scl(:,1) morph_scl(:,2)],sLabels,'k',12);
% 
% Directions of variation:
figure;
set(gcf,'color','w'); % set bg color to white
f_biplotPca2(morph.evects,1,0,mLabels);
f_biplotPca2(perform.evects(:,1:2)*H,1,0,pLabels);
box on;
