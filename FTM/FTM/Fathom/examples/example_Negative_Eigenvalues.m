% Example of handling negative eigenvalues
% 
% by David L. Jones, Mar-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% Jun-2012: updated to reflect changes in f_pcoa 

% File: '.../examples/legendre_113.mat'

% -----NOTES:-----
% This file loosely reproduces the RDA analysis of coral reef fish abundances in
% Table 11.3 of Legendre & Legendre (1998), but employs a Bray-Curtis
% distance-based approach and focuses only on the qualitative explanatory
% variable 'substrate type' and not depth.

% -----Variables:-----
% y     = abundance of 6 species of fishes from 10 sites
% y_txt = corresponding column labels
% x     = column vector of integers specifying substrate type associated with
%        each site (1=coral, 2=sand, 3=other)
% x_txt = cell array of corresponding category labels
% depth = water depth associated with each site

% -----ANALYSIS:-----
% Load data file:
load legendre_113.mat

% Square-root transform the data to reduce the effect of more abundance species:
y_2 = f_normal(y,'2');

% Bray-Curtis dissimilarities among samples:
disBC = f_dis(y_2,'bc');

% Dummy code the substrate type (trim last column to avoid a singular matrix):
x_dum = f_dummy(x,1);

% PCoA:
U_0 = f_pcoa(disBC,0,1,0); % discard negative eigenvalues
U_1 = f_pcoa(disBC,0,1,1); % keep negative eigenvalues
U_2 = f_pcoa(disBC,0,1,2); % correct negative eigenvalues

% db-RDA (discard neg eigenvalues):
rda_0 = f_rda(U_0.scores,x_dum,0,1000,1);
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 4.0863    p    =  0.00100 
% R2 = 0.5386   R2adj =  0.40683 
% No. of permutations = 1000 
% --------------------------------------------------

% db-RDA (keep neg eigenvalues):
rda_1 = f_rda(U_1.scores,x_dum,0,1000,1);
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 3.9770    p    =  0.00100 
% R2 = 0.5319   R2adj =  0.39815 
% No. of permutations = 1000 
% --------------------------------------------------

% db-RDA (correct neg eigenvalues):
rda_2 = f_rda(U_2.scores,x_dum,0,1000,1);
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 3.7985    p    =  0.00100 
% R2 = 0.5205   R2adj =  0.38344 
% No. of permutations = 1000 
% --------------------------------------------------

% Compare with NP-MANOVA:
result = f_npManova(disBC,x,1000,1);
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'         'F'         'p'    
%     'factor 1'    [ 2]    [1.2553]    [0.62766]    [4.2018]    [0.001]
%     'residual'    [ 7]    [1.0457]    [0.14938]    [   NaN]    [  NaN]
%     'total'       [ 9]    [ 2.301]    [    NaN]    [   NaN]    [  NaN]
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 


% -> NP-MANOVA provides a true measure of the SS Total, whereas in db-RDA with a
% semi-metric dissimilarity metric the SS Total is inflated as the SS associated
% with the negative eigenvalues needs to be subtracted from the SS associated
% with the positive eigenvalues


% Use a better method of db-RDA which directly (& correctly) estimates SS's:
result_DB = f_rdaDB(disBC,size(y_2,2),x_dum,0,1000,1);
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 4.2018    p    =  0.00100 
% R2 = 0.5456   R2adj =  0.41572 
% No. of permutations = 1000 
% --------------------------------------------------

% RDA distance biplot with vectors for X,Y:
result_DB.X = f_dummy(x,0); % replace trimmed X with untrimmed version
f_rdaPlot(result_DB,y_2,0,[1 1/15],0,y_txt,x_txt);

% Customize plot:
axisVar = axis;                 % get current axes bounds
axis([-0.5 0.75 axisVar(3:4)]); % adjust x-axis
