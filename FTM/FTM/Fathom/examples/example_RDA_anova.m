% Example of using RDA to perform a 2-way ANOVA
% 
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% updated: Mar-2010

% File: '.../examples/rats.mat'

% -----NOTES:-----
% Data is from Table 11.1 of Sokal & Rohlf (1995) and is also available as
% 'sr_11p1.mat' in FATHOM's data folder

% Load data file:
load rats.mat

% Perform ANOVA using RDA:
anova = f_rdaAnova(food,[sex fat],1000,1);
% > 22
% 
% ==================================================
%     Permutation-based (M)ANOVA via RDA:
% --------------------------------------------------
%     'Source'        'df'    'SS'        'MS'        'F'         'p'    
%     'factor 1'      [ 1]    [3780.8]    [3780.8]    [2.5925]    [0.182]
%     'factor 2'      [ 1]    [ 61204]    [ 61204]    [41.969]    [0.001]
%     'factor 1x2'    [ 1]    [918.75]    [918.75]    [  0.63]    [0.403]
%     'residual'      [ 8]    [ 11667]    [1458.3]    [   NaN]    [  NaN]
%     'total'         [11]    [ 77570]    [   NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 


% Perform same analysis using NP-MANOVA:
dis    = f_dis(food,'euc');
result = f_npManova(dis,[sex fat],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% 
% Select model...[O will cancel]
% 21
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'        'MS'        'F'         'p'    
%     'factor 1'      [ 1]    [3780.8]    [3780.8]    [2.5925]    [0.141]
%     'factor 2'      [ 1]    [ 61204]    [ 61204]    [41.969]    [0.001]
%     'factor 1x2'    [ 1]    [918.75]    [918.75]    [  0.63]    [0.426]
%     'residual'      [ 8]    [ 11667]    [1458.3]    [   NaN]    [  NaN]
%     'total'         [11]    [ 77570]    [   NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication)
