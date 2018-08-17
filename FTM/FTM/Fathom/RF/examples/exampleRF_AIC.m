% Example Random Forest Variable Selection
% 
% -----Author:-----
% by David L. Jones, Sep-2009
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% File: '.../examples/iris_R.mat'

% Load Fisher's Iris data:
load iris_R.mat X Y

% Variable selection via stepwise forward addition:
% [best,cond] = f_RFaic(X,Y,nTree,BIC,verb,cut,xLabels);
best = f_RFaic(X,Y,[],0,1);
% ==================================================
% AIC-based stepwise forward selection (RANDOM FOREST)
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'ERR'         'AIC'        'wts'       'delta'          'ratio'         'var' 
%     [0.046667]    [-455.63]    [     0]    [          1]    [         1]    '4'   
%     [    0.08]    [-374.78]    [80.849]    [ 2.7782e-18]    [3.5995e+17]    '3'   
%     [ 0.29333]    [-179.89]    [275.74]    [ 1.3286e-60]    [7.5266e+59]    '1'   
%     [     0.5]    [ -99.89]    [355.74]    [ 5.6592e-78]    [ 1.767e+77]    '2'   
%     [       1]    [  2.027]    [457.65]    [4.1847e-100]    [2.3897e+99]    'none'
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'ERR'         'AIC'        'wts'    'deltaN'    'var'     'idx'
%     [0.046667]    [-455.63]    [  1]    [457.65]    '4'       [  4]
%     [0.033333]    [-504.02]    [  1]    [48.388]    '3'       [  3]
%     [0.033333]    [-504.02]    [  1]    [     0]    'none'       []
% --------------------------------------------------
% 
% # Trees = 1000
% ERR     = classification error 
% AIC     = corrected AIC 
% deltaN  = delta associated with NO variable addition 
% wts     = AIC weights 
% var     = variable labels 
% idx     = index to selected variables 
% 
% (Note: ERR in Marginal Tests are CUMULATIVE) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            GENERATE RANDOM FOREST:            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rf = f_RFclass(X(:,best.idx),Y,[],[],0,0,'raw',1);
% =======================================================
% =====     INPUT PARAMETERS FOR RANDOM FOREST:     =====
% =======================================================
% nTree          = 1000
% mTry           = 1
% X: # obs = 150, # var's = 2 
% # Y classes    = 3 
% addclass       = 0
% size(ncat)     = 1
% size(ncat)     = 2
% maxcat         = 1
% size(sampsize) = 1
% size(sampsize) = 1
% sampsize[0]    = 150
% stratify       = 0
% imp            = 0
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
% Ave. error rate (after 1000 iterations) = 0.0333
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    1           100.00 % 
%    2            94.00 % 
%    3            96.00 % 
% 
% 
% Total Correct  = 96.67 % 
% Total Error    = 3.33 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      1      2      3 
%      1 100.00   0.00   0.00 
%      2   0.00  94.00   6.00 
%      3   0.00   4.00  96.00 
% ==================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   COMPARE WITH RDA:   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
f_rdaAIC(f_designMatrix(Y,0),X,0,1);
% ==================================================
% AIC-based stepwise forward selection (RDA)         
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'       'R2'         'R2adj'      'AIC'        'wts'       'delta'         'ratio'         'var' 
%     [52.931]    [0.47069]    [0.46711]    [-152.16]    [     0]    [   0.70674]    [         1]    '3'   
%     [53.556]    [0.46444]    [0.46082]    [ -150.4]    [1.7592]    [   0.29326]    [      2.41]    '4'   
%     [69.065]    [0.30935]    [0.30469]    [-112.26]    [39.907]    [ 1.526e-09]    [4.6313e+08]    '1'   
%     [79.961]    [0.20039]    [0.19499]    [-90.283]    [61.881]    [ 2.582e-14]    [2.7371e+13]    '2'   
%     [   100]    [    NaN]    [    NaN]    [-58.793]    [93.371]    [3.7491e-21]    [1.8851e+20]    'none'
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'         'R2adj'      'AIC'        'wts'        'deltaN'    'var'     'idx'
%     [52.931]    [0.47069]    [0.46711]    [-152.16]    [0.70674]    [93.371]    '3'       [  3]
%     [44.005]    [0.55995]    [0.55397]    [-177.79]    [0.99753]    [25.623]    '2'       [  2]
%     [40.504]    [0.59496]    [0.58663]    [-188.11]    [0.98852]    [10.321]    '4'       [  4]
%     [40.504]    [    NaN]    [    NaN]    [-188.11]    [0.70815]    [     0]    'none'       []
% --------------------------------------------------
% 
% RSS   = residual sum-of-squares 
% R2    = fraction of total variance explained 
% R2adj = fraction of adjusted total variance explained 
% AIC   = corrected AIC 
% deltaN = delta associated with NO variable addition 
% wts   = AIC weights 
% var   = variable labels 
% idx   = index to selected variables 
% 
% (Note: RSS, R2, and R2adj in Marginal tests are CUMULATIVE) 


% -> note there are differences in the variables selected by the two methods,
%    principally because RDA optimized fit of the model, while RF optimized
%    classification success





