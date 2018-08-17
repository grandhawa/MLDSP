% Example Random Forest: Small Sample Size
% 
% -----Author:-----
% by David L. Jones, Nov-2010
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     FISHER'S IRIS DATA:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Fisher's Iris data:
load iris_R.mat

% Extract a small random subset:
[sX sY] = f_grpResample(X,Y,5,0);

% Bootstrap small sample:
B = f_grpBoot(sX,sY,0);

% Generate Random Forest:
% USAGE: model = f_RFclass(X,Y,nTree,mTry,imp,sim,'stnd',verb,X_txt,opt);
model = f_RFclass(B,sY,[],[],1,0,'raw',1);

% Generating 1000 trees with 2 random subsets of variables
% for 3 classes...
% 
% =======================================================
% =====     INPUT PARAMETERS FOR RANDOM FOREST:     =====
% =======================================================
% nTree          = 1000
% mTry           = 2
% X: # obs = 15, # var's = 4 
% # Y classes    = 3 
% addclass       = 0
% size(ncat)     = 1
% size(ncat)     = 4
% maxcat         = 1
% size(sampsize) = 1
% size(sampsize) = 1
% sampsize[0]    = 15
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
%    1           100.00 % 
%    2           100.00 % 
%    3           100.00 % 
% 
% 
% Total Correct  = 100.00 % 
% Total Error    = 0.00 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      1      2      3 
%      1 100.00   0.00   0.00 
%      2   0.00 100.00   0.00 
%      3   0.00   0.00 100.00 
% ==================================================


% Calculate classification success using TEST data:
Y_hat = f_RFclassPredict(X,model,Y,0,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        SATIMAGE DATA:       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data
load satimage.mat

% Extract a small random subset:
[sX_train sY_train] = f_grpResample(X_train,Y_train,10,0);

% Generate Random Forest using TRAINING data:
[model,err] = f_RFclass(sX_train,sY_train,[],[6],0,0,'stnd',1);
% 
% Ave. error rate (after 1000 iterations) = 0.2333
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    1            60.00 % 
%    2           100.00 % 
%    3            80.00 % 
%    4            80.00 % 
%    5            60.00 % 
%    6            80.00 % 
% 
% 
% Total Correct  = 76.67 % 
% Total Error    = 23.33 % 
% Prior prob     = class size 

% Calculate classification success using TEST data:
Y_hat = f_RFclassPredict(X_test,model,Y_test,0,1);
% % 
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%     Classification Success Using TEST Data: 
% --------------------------------------------------
% Group        Corrrect  
%    1            80.48 % 
%    2            88.39 % 
%    3            85.39 % 
%    4            76.78 % 
%    5            81.01 % 
%    6            67.23 % 
% 
% 
% Total Correct  = 78.90 % 
% Total Error    =  21.10 % 
% Prior prob     = group size 



% Bootstrap a larger sample:
[B bGrp] = f_grpResample(sX_train,sY_train,100,1,0);

% Generate Random Forest using TRAINING data:
[modelB,errB] = f_RFclass(B,bGrp,[],[6],0,1,'stnd',1);
% Ave. error rate (after 1000 iterations) = 0.0000
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    1           100.00 % 
%    2           100.00 % 
%    3           100.00 % 
%    4           100.00 % 
%    5           100.00 % 
%    6           100.00 % 
% 
% 
% Total Correct  = 100.00 % 
% Total Error    = 0.00 % 
% Prior prob     = class size 

% Calculate classification success using TEST data:
Y_hat = f_RFclassPredict(X_test,modelB,Y_test,0,1);
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%     Classification Success Using TEST Data: 
% --------------------------------------------------
% Group        Corrrect  
%    1            78.74 % 
%    2            87.05 % 
%    3            82.87 % 
%    4            76.30 % 
%    5            83.97 % 
%    6            57.23 % 
% 
% 
% Total Correct  = 75.80 % 
% Total Error    =  24.20 % 
% Prior prob     = group size 
