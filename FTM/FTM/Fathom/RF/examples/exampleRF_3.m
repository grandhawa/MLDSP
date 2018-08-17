% Example Random Forest 3
% 
% -----Author:-----
% by David L. Jones, Sep-2009
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----References:-----
% http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#micro

% Load the microarray lymphoma data set (81 observations, 4682 variables
% corresponding to gene expressions, and 3 classes):
load microarray.mat

% Create Random Forest:
[model,err] = f_RFclass(X,Y,[1000],[150],1,0,'raw',1);
% 
% Ave. error rate (after 1000 iterations) = 0.0123
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    0            96.55 % 
%    1           100.00 % 
%    2           100.00 % 
% 
% 
% Total Correct  = 98.77 % 
% Total Error    =  1.23 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      0      1      2 
%      0  96.55   3.45   0.00 
%      1   0.00 100.00   0.00 
%      2   0.00   0.00 100.00 
% 
% ==================================================


% List the 25 most important variables, sorted by their z-scores:
idx = f_RFimp(model,25);
% 
% ==================================================
% =====      Variable Importance Measures:     =====
% ==================================================
% col 1: variable label
% col 2: meanAcc
% col 3: zScore
% col 2: zP
% ---------------------
%  689   1.666   1.121   0.131
%  666   1.359   1.043   0.149
%  668   1.281   0.998   0.159
%  663   0.894   0.894   0.186
%  667   1.086   0.885   0.188
%  682   0.798   0.846   0.199
% 3622   0.666   0.814   0.208
%  685   0.724   0.737   0.231
%  650   0.555   0.723   0.235
% 1080   0.513   0.715   0.237
%  686   0.654   0.701   0.242
%  897   0.608   0.699   0.242
%  697   0.549   0.697   0.243
%  669   0.410   0.663   0.254
% 3631   0.453   0.650   0.258
% 3405   0.352   0.650   0.258
% 1427   0.457   0.634   0.263
%  895   0.426   0.627   0.265
%  723   0.434   0.618   0.268
%  679   0.415   0.617   0.268
% 3531   0.360   0.612   0.270
%  656   0.305   0.592   0.277
%  681   0.434   0.592   0.277
% 1055   0.350   0.590   0.278
% 1104   0.407   0.589   0.278
% ==================================================



% Create a new Random Forest based on top 15 variables (include labels):
idx = idx(1:15); % index to top 15
[model_2,err_2] = f_RFclass(X(:,idx),Y,[1000],[],1,0,'raw',1);
% 
% Generating 1000 trees with 3 random subsets of variables
% for 3 classes...
% 
% =======================================================
% =====     INPUT PARAMETERS FOR RANDOM FOREST:     =====
% =======================================================
% nTree          = 1000
% mTry           = 3
% X: # obs = 81, # var's = 15 
% # Y classes    = 3 
% addclass       = 0
% size(ncat)     = 1
% size(ncat)     = 15
% maxcat         = 1
% size(sampsize) = 1
% size(sampsize) = 1
% sampsize[0]    = 81
% stratify       = 0
% imp            = 1
% sim            = 0
% oob_prox       = 0
% strata         = 1
% ipi            = 0
% classwt        = 1.0000
% classwt        = 1.0000
% classwt        = 1.0000
% cutoff         = 0.3333
% cutoff         = 0.3333
% cutoff         = 0.3333
% nodesize       = 1.0000
% =======================================================
% (Using raw data for X) 
% Ave. error rate (after 1000 iterations) = 0.0000
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    0           100.00 % 
%    2           100.00 % 
%    1           100.00 % 
% 
% 
% Total Correct  = 100.00 % 
% Total Error    = 0.00 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      0      2      1 
%      0 100.00   0.00   0.00 
%      2   0.00 100.00   0.00 
%      1   0.00   0.00 100.00 
% ==================================================

% -> Note the classification error using only 15 variables is better than one
%    using the entire suite of 4682 variables.

% List the importance measures for all 15 variables used:
f_RFimp(model_2);
% 
% ==================================================
% =====      Variable Importance Measures:     =====
% ==================================================
% col 1: variable label
% col 2: meanAcc
% col 3: SD
% col 4: zP
% ---------------------
%  4   0.067   0.030   0.013
% 14   0.064   0.024   0.004
% 13   0.063   0.023   0.003
%  3   0.063   0.029   0.014
%  2   0.063   0.029   0.016
%  6   0.060   0.023   0.005
%  1   0.058   0.029   0.021
%  7   0.053   0.023   0.010
%  5   0.044   0.024   0.030
%  9   0.044   0.020   0.013
%  8   0.024   0.018   0.093
% 11   0.019   0.015   0.108
% 15   0.018   0.014   0.103
% 12   0.017   0.014   0.118
% 10   0.017   0.014   0.120
% ==================================================
