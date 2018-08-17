% Examples of performing MANOVA via db-RDA
% 
% by David L. Jones, Oct-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 

% 1) Example of a One-Way Model I ANOVA with replication (balanced design). This
% is data from Box 9.5 of Sokal & Rohlf (1995) and has 1 response variable
% (age), and 1 fixed factor (clone).
% 
% Load the data file:
load sr_09p5.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
yDis      = f_dis(age,'euc');
result_01 = f_rdaDB_manova(yDis,clone,1000,1);
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'           'MS'           'F'           'p'    
%     'factor 1'    [ 1]    [0.0064286]    [0.0064286]    [0.014063]    [0.832]
%     'residual'    [12]    [   5.4857]    [  0.45714]    [     NaN]    [  NaN]
%     'total'       [13]    [   5.4921]    [      NaN]    [     NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 


% 2) Example of a Two-way Model I ANOVA with replication. This data is from
% Table 11.1 of Sokal & Rohlf (1995) and has 1 response variable (food) and 2
% fixed factors: fat and sex.
% 
% Load the data file:
load sr_11p1.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
yDis      = f_dis(food,'euc');
result_02 = f_rdaDB_manova(yDis,[sex fat],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% Select model...[O will cancel]
% 21
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'        'MS'        'F'         'p'    
%     'factor 1'      [ 1]    [3780.8]    [3780.8]    [2.5925]    [ 0.18]
%     'factor 2'      [ 1]    [ 61204]    [ 61204]    [41.969]    [0.001]
%     'factor 1x2'    [ 1]    [918.75]    [918.75]    [  0.63]    [0.471]
%     'residual'      [ 8]    [ 11667]    [1458.3]    [   NaN]    [  NaN]
%     'total'         [11]    [ 77570]    [   NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication)
