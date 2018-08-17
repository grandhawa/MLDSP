% Example of using Random Forest for Maximum Likelihood Estimation
% 
% -----Author:-----
% by David L. Jones, May-2013
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/hisea.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load HISEA data:
load hisea.mat

% Create group labels:
grp_txt = repmat({'ak1202b'},size(X.stock,1),1);
idx2 = find(X.stock==2);
idx3 = find(X.stock==3);
grp_txt(idx2) = {'nas1202'};
grp_txt(idx3) = {'ske1202'};
clear idx2 idx3;

% Create a Random Forest using the baseline data:
rf = f_RFclass(X.dat,X.stock,[],[],0,0,'raw',1,[]);
% 
% Generating 1000 trees with 2 random subsets of variables
% for 3 classes...
% 
% =======================================================
% =====     INPUT PARAMETERS FOR RANDOM FOREST:     =====
% =======================================================
% nTree          = 1000
% mTry           = 2
% X: # obs = 598, # var's = 4 
% # Y classes    = 3 
% addclass       = 0
% size(ncat)     = 1
% size(ncat)     = 4
% maxcat         = 1
% size(sampsize) = 1
% size(sampsize) = 1
% sampsize[0]    = 598
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
% Ave. error rate (after 1000 iterations) = 0.2375
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    1            73.00 % 
%    2            71.86 % 
%    3            83.92 % 
% 
% 
% Total Correct  = 76.25 % 
% Total Error    = 23.75 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      1      2      3 
%      1  73.00  22.00   5.00 
%      2  15.58  71.86  12.56 
%      3   4.52  11.56  83.92 
% ==================================================

% Estimate mixture proportions using a Maximum Likelihood approach:
mle = f_RFmle(Y.dat,rf,grp_txt,1);
% 
% ==================================================
%      RF - MAXIMUM LIKELIHOOD ESTIMATION
%          Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.2863 :ak1202b 
% Group 2 = 0.3209 :nas1202 
% Group 3 = 0.3929 :ske1202 
% 
% # iterations for EM convergence = 5 
% ==================================================

% Determine the agreement in classifications determined by the MLE approach vs.
% that based on majority rules voting by the Random Forest (in %):
sum(mle.grpM == mle.Y_hat)/size(mle.grpM,1) * 100
% ans = 98.462

% Show the overall mean level of confidence of the classification by
% Random Forest:
mean(mle.marg)
% ans = 0.57903



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/codmix.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data from the mixFish package for R:
load codmix.mat;
% 
% Variables:
% raw = structure of data with the following fields:
%  .dat        = otolith microchemical + genetic data
%  .txt        = corresponding variable labels
%  .season     = 1 (S96) or 2 (J96)
%  .season_txt = corresponding text labels
%  .stock      = 1 ('3Pn4RS'), 2('3Ps'), 3 ('4T'), 4 ('4Vn'), 5 ('4Vs')
%  .stock_txt  = corresponding text labels              

% Parse only the otolith data:
idxX  = find(raw.season==1); % Training set
X.dat = raw.dat(idxX,1:5);
X.grp = raw.stock(idxX);
X.txt = raw.stock_txt(idxX);
idxY  = find(raw.season==2); % Sample data
Y.dat = raw.dat(idxY,1:5);
% Clean up:
clear idxX idxY;

% Create a Random Forest using the baseline data:
rf = f_RFclass(X.dat,X.grp,[],[],0,0,'raw',1,[]);
% 
% 
% Generating 1000 trees with 2 random subsets of variables
% for 4 classes...
% 
% =======================================================
% =====     INPUT PARAMETERS FOR RANDOM FOREST:     =====
% =======================================================
% nTree          = 1000
% mTry           = 2
% X: # obs = 200, # var's = 5 
% # Y classes    = 4 
% addclass       = 0
% size(ncat)     = 1
% size(ncat)     = 5
% maxcat         = 1
% size(sampsize) = 1
% size(sampsize) = 1
% sampsize[0]    = 200
% stratify       = 0
% imp            = 0
% sim            = 0
% oob_prox       = 0
% strata         = 1
% ipi            = 0
% classwt        = 1.0000
% classwt        = 1.0000
% classwt        = 1.0000
% classwt        = 1.0000
% cutoff         = 0.2500
% cutoff         = 0.2500
% cutoff         = 0.2500
% cutoff         = 0.2500
% nodesize       = 1.0000
% =======================================================
% (Using raw data for X) 
% Ave. error rate (after 1000 iterations) = 0.4900
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    3            57.89 % 
%    1            34.48 % 
%    5            63.04 % 
%    2            45.00 % 
% 
% 
% Total Correct  = 51.00 % 
% Total Error    = 49.00 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      3      1      5      2 
%      3  57.89  30.26  10.53   1.32 
%      1  53.45  34.48   6.90   5.17 
%      5  19.57   8.70  63.04   8.70 
%      2   0.00  20.00  35.00  45.00 
% ==================================================

% Estimate mixture proportions using a RF-MLE approach:
mle = f_RFmle(Y.dat,rf,X.txt,1);
% 
% ==================================================
%      RF - MAXIMUM LIKELIHOOD ESTIMATION
%          Mixing Proportions:
% --------------------------------------------------
% Group 3 = 0.4247 :4T 
% Group 1 = 0.0000 :3Pn4RS 
% Group 5 = 0.5455 :4Vs 
% Group 2 = 0.0298 :3Ps 
% 
% # iterations for EM convergence = 16 
% ==================================================

% -> Since the classification success rate of this RF performs relatively poorly,
%    with only a 51% overall success rate, it might be better to use a
%    different method to estimate mixture proportions


%%%%%%%%%%%%%%%%%%%%%%%%
%     CAP - MLE:       %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Determine the optimal value of 'm' for a CAP model:
f_capOptimal(X.dat,'euc',X.grp,1,1)
% ==================================================
%  Diagnostics for CAP - db-CDA:
% --------------------------------------------------
% m:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:
% 1  99.9993  3.5454 0.3288 0.0000 0.0000   27.00% <- optimal value of m = 1
% 2  99.9998  3.4818 0.3337 0.1293 0.0000   19.50% 
% 3  99.9999  3.5466 0.3905 0.1980 0.0025   22.00% 
% 4  100.0000  3.5603 0.3905 0.2000 0.0038   22.00% 
% 5  100.0000  3.5970 0.3908 0.2458 0.0047   20.50% 
% --------------------------------------------------
% m       = # PCoA axes retained
% propG   = proportion of yDis explained by m PCoA axes
% RSS     = leave-one-out residual sums-of-squares
% d_1^2   = squared canonical correlation for axis 1
% Correct = leave-one-out classification success
% 
% Central Tendency = spatial median 
% Optimal value of m for CDA may be 1
% ==================================================
% 
% 
% Estimate mixture proportions using a CAP-MLE approach:
cap_mle = f_capMLE(X.dat,'euc',X.grp,X.txt,Y.dat,1,1,1);
% 
% ==================================================
%      CAP - MAXIMUM LIKELIHOOD ESTIMATION
%          Mixing Proportions:
% --------------------------------------------------
% Group 3 = 0.9170 :4T 
% Group 1 = 0.0000 :3Pn4RS 
% Group 5 = 0.0823 :4Vs 
% Group 2 = 0.0007 :3Ps 
% 
% # iterations for EM convergence = 37 
% ==================================================
% 
% -> the CAP model was a poorly performing classifier, so it might not be
%    the best method to estimate mixture proportions.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/iris.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Fisher's Iris data:
load iris.mat

% Partition the data into separate BASELINE and MIXED sets
idx = f_shuffle(1:size(grps,1));
B      = iris(idx(1:125),:);
B_grps = grps(idx(1:125));
M      = iris(idx(126:end),:);
M_grps = grps(idx(126:end));

% Sort the BASELINE data (required by CAP):
[nul,idxS] = sortrows(B_grps);
B          = B(idxS,:);
B_grps     = B_grps(idxS);
clear idx idxS nul;

% Save the data:
fname = 'iris_test';
saver;



%%%%%%%%%%%%%%%%%%%%%%%%
%      RF - MLE:       %
%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Random Forest using the baseline data:
rf = f_RFclass(B,B_grps,[],[],0,0,'raw',1,[]);
% 
% Generating 1000 trees with 2 random subsets of variables
% for 3 classes...
% 
% =======================================================
% =====     INPUT PARAMETERS FOR RANDOM FOREST:     =====
% =======================================================
% nTree          = 1000
% mTry           = 2
% X: # obs = 125, # var's = 4 
% # Y classes    = 3 
% addclass       = 0
% size(ncat)     = 1
% size(ncat)     = 4
% maxcat         = 1
% size(sampsize) = 1
% size(sampsize) = 1
% sampsize[0]    = 125
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
% Ave. error rate (after 1000 iterations) = 0.0480
% ==================================================
%        RANDOM FOREST 'Majority Rules' 
%   Internal Cross-validation Classification Success: 
% --------------------------------------------------
% Class        Corrrect  
%    1           100.00 % 
%    2            95.35 % 
%    3            90.70 % 
% 
% 
% Total Correct  = 95.20 % 
% Total Error    = 4.80 % 
% Prior prob     = class size 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% class:      1      2      3 
%      1 100.00   0.00   0.00 
%      2   0.00  95.35   4.65 
%      3   0.00   9.30  90.70 
% ==================================================

% Estimate mixture proportions using a RF-MLE approach:
mle = f_RFmle(M,rf,[],1);
% 
% ==================================================
%      RF - MAXIMUM LIKELIHOOD ESTIMATION
%          Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.4426 
% Group 2 = 0.2387 
% Group 3 = 0.3187 
% 
% # iterations for EM convergence = 3 
% ==================================================

% Compare to the true mixture proportions:
for i=1:3
   mProp(i) = sum(M_grps==i)/numel(M_grps);
end
mProp
% mProp = 0.44 0.28 0.28




%%%%%%%%%%%%%%%%%%%%%%%%
%     CAP - MLE:       %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Determine the optimal value of 'm' for a CAP model:
f_capOptimal(B,'euc',B_grps,1,1);
% 
% ==================================================
%  Diagnostics for CAP - db-CDA:
% --------------------------------------------------
% m:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:
% 1  92.2699  1.8798 0.9274 0.0000 0.0000   93.60% <- optimal value of m
% 2  97.6720  1.7207 0.9620 0.1319 0.0000   92.80% 
% 3  99.4436  1.8261 0.9666 0.2472 0.0000   92.80% 
% 4  100.0000  1.8434 0.9678 0.2591 0.0000   93.60% 
% --------------------------------------------------
% m       = # PCoA axes retained
% propG   = proportion of yDis explained by m PCoA axes
% RSS     = leave-one-out residual sums-of-squares
% d_1^2   = squared canonical correlation for axis 1
% Correct = leave-one-out classification success
% 
% Central Tendency = spatial median 
% Optimal value of m for CDA may be 1
% ==================================================
% 
% Estimate mixture proportions using a CAP-MLE approach:
cap_mle = f_capMLE(B,'euc',B_grps,[],M,1,[1],1);

% ==================================================
%      CAP - MAXIMUM LIKELIHOOD ESTIMATION
%          Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.0003 
% Group 2 = 0.9997 
% Group 3 = 0.0000 
% 
% # iterations for EM convergence = 39 
% ==================================================


% cap_mle = 
% 
%     theta: [0.00028431 0.99972 9.8072e-11]
%      grpM: [25x1 double]
%        PP: [25x3 double]
%       grp: [1 2 3]
%       txt: {'1'  '2'  '3'}
%     nIter: 39
%       raw: [0.44 0.32 0.24]
%     PPcap: [25x3 double]
%      marg: [25x1 double]
%      grpW: [25x1 double]

% -> theta is WAY off, but the raw proportions are much closer



% MLE of mixture proportions (UNEQUAL COV):
resultIris = f_mle(iris,iris,grps,[],1,2,100,250);
% 
% ==================================================
%         MAXIMUM LIKELIHOOD ESTIMATION
%             Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.3333 (0.0385) 
% Group 2 = 0.3254 (0.0409) 
% Group 3 = 0.3413 (0.0432) 
% 
% PDF = MVN + unequal covariances
% 
% # iterations for EM convergence = 2 
% (Bootstrapped SD's are in parentheses)
% 
% --------------------------------------------------
% Bootstrapped correlations among Groups:
% Group 1 vs. 2 = -0.4095 
% Group 1 vs. 3 = -0.5033 
% Group 2 vs. 3 = -0.5823 
% 
% # iterations for Bootstrap = 250 
% 
% NOTE: Strong neg. correlations between groups suggests
% the MLE procedure considers these groups similar
% ==================================================
% 
% resultIris = 
%     theta: [0.33333 0.32541 0.34126]
%        SD: [0.038524 0.040948 0.04323]
%       cor: [3x3 double]
%      grpM: [150x1 double]
%        PP: [150x3 double]
%       grp: [1 2 3]
%       txt: {'1'  '2'  '3'}
%     nIter: 2
%       pdf: 'MVN + unequal covariances'
%       raw: [0.33333 0.32667 0.34]
