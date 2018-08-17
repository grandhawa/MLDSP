% Example for Canonical Discriminant Analysis
% 
% by David L. Jones, May-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% Jan-2011: updated
% Jun-2012: added f_chanceClass

% File: '.../examples/iris.mat'

% Load Fisher's Iris data:
load iris.mat

% Perform CDA:
result = f_cda(iris,grps,2,1000,1);
% 
% Permuting the data 999 times...
% ==================================================
% CANONICAL DISCRIMINANT ANALYSIS:
%  Discriminant Functions (method = 2) 
% --------------------------------------------------
% Trace stat          = 2387.0830  p =  0.00100 
% Greatest root stat  = 2366.1068  p =  0.00100 
% No. of permutations = 1000 
% Wilks' lambda       = 0.0234 
% --------------------------------------------------
% Canonical Eigenvalues:
%   2366.1068  20.9762
% Fraction of AMONG-GROUP variance explained:
% ------------------------------
% Canonical axes: 
%   0.9912  0.0088
%   cumulative:
%   0.9912  1.0000
% ==================================================

% Make CDA Plot:
f_cdaPlot(result,0,0,0,0,(1),[0.04],{'one','two','three','four'});

% Perform 'LOO' cross-validation:
[err,PP] = f_cdaCV(iris,grps);
% ==================================================
%             F_CDA CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Corrrect  
%    1            100.0 % 
%    2             96.0 % 
%    3             98.0 % 
% 
% 
% Total Correct  = 98.00 % 
% Total Error    =  2.00 % 
% Class. method  = centroid 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      3      
%      1  100.0    0.0    0.0 
%      2    0.0   96.0    4.0 
%      3    0.0    2.0   98.0 
% 
% ==================================================

% Test the significance of the observed classification success rate:
f_chanceClass(grps,1,1000,1);
% 
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1            11.11 % 
%    2            11.11 % 
%    3            11.11 % 
% 
% Total Correct = 33.33 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 33.40 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% -> Reject the null hypothesis that the observed re-classification success
%    rate is no better than that expected by chance alone, at alpha = 0.05.
%    