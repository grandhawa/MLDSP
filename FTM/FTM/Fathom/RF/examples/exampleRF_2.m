% Example Random Forest
% 
% -----Author:-----
% by David L. Jones, Aug-2009
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% File: '.../examples/satimage.mat'

% -----References:-----
% http://www.stat.berkeley.edu/~breiman/RandomForests/cc_manual.htm#satimage

% The 'satimage' data consists of a training set with 4435 observations, 36
% variables, and 6 classes with a corresponding test set of 2000 observations.

% Load the data
load satimage.mat

% Generate Random Forest using TRAINING data:
[model,err] = f_RFclass(X_train,Y_train,[50],[6],0,0,'raw',1);

% Ave. error rate (after 50 iterations) = 0.0954
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    1            97.29 % 
%    2            96.87 % 
%    3            95.53 % 
%    4            58.80 % 
%    5            85.74 % 
%    6            90.56 % 
% 
% 
% Total Correct  = 90.46 % 
% Total Error    = 9.54 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      1      2      3      4      5      6 
%      1  97.29   0.19   1.49   0.00   1.03   0.00 
%      2   0.00  96.87   0.21   0.63   1.25   1.04 
%      3   0.94   0.10  95.53   2.39   0.00   1.04 
%      4   0.96   0.48  21.20  58.80   0.72  17.83 
%      5   4.47   1.06   0.21   1.28  85.74   7.23 
%      6   0.00   0.00   1.73   5.68   2.02  90.56 
% ==================================================


% Calculate classification success using TEST data:
Y_hat = f_RFclassPredict(X_test,model,Y_test,0,1);

% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%     Classification Success Using TEST Data: 
% --------------------------------------------------
% Group        Corrrect  
%    1            99.13 % 
%    2            96.43 % 
%    3            94.46 % 
%    4            63.51 % 
%    5            89.03 % 
%    6            89.79 % 
% 
% 
% Total Correct  = 90.75 % 
% Total Error    =  9.25 % 
% Prior prob     = group size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group:      1      2      3      4      5      6 
%      1  99.13   0.00   0.22   0.00   0.65   0.00 
%      2   0.00  96.43   0.89   0.45   1.34   0.89 
%      3   0.76   0.00  94.46   3.02   0.25   1.51 
%      4   0.00   0.47  14.22  63.51   0.95  20.85 
%      5   2.95   1.27   0.00   0.42  89.03   6.33 
%      6   0.00   0.00   2.77   5.53   1.91  89.79 
% 
% ==================================================












