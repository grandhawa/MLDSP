% Example for Maximum Likelihood Estimation
% 
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/hisea.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load HISEA data:
load hisea.mat

% MLE of mixture proportions (MVN + EQUAL COV):
% (A) No bootstrap estimation of SD or correlation:
% 
result = f_mle(Y.dat,X.dat,X.stock,[],1,1)
% 
% ==================================================
%         MAXIMUM LIKELIHOOD ESTIMATION
%             Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.2424 (NaN) 
% Group 2 = 0.3702 (NaN) 
% Group 3 = 0.3874 (NaN) 
% 
% PDF = MVN + equal covariance
% 
% # iterations for EM convergence = 4 
% (Bootstrapped SD's are in parentheses)
% 
% ==================================================
% 
% result = 
%     theta: [0.24235 0.37024 0.3874]
%        SD: [NaN NaN NaN]
%       cor: NaN
%      grpM: [195x1 double]
%        PP: [195x3 double]
%       grp: [1 2 3]
%       txt: {'1'  '2'  '3'}
%     nIter: 4
%       pdf: 'MVN + equal covariance'
%       raw: [0.26154 0.35385 0.38462]
 
 
% (B) With bootstrap estimation of SD & correlation:
% 
result = f_mle(Y.dat,X.dat,X.stock,[],1,1,100,500)
% 
% ==================================================
%         MAXIMUM LIKELIHOOD ESTIMATION
%             Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.2424 (0.0405) 
% Group 2 = 0.3702 (0.0485) 
% Group 3 = 0.3874 (0.0440) 
% 
% PDF = MVN + equal covariance
% 
% # iterations for EM convergence = 4 
% (Bootstrapped SD's are in parentheses)
% 
% --------------------------------------------------
% Bootstrapped correlations among Groups:
% Group 1 vs. 2 = -0.5233 
% Group 1 vs. 3 = -0.3437 
% Group 2 vs. 3 = -0.6204 
% 
% # iterations for Bootstrap = 500 
% 
% NOTE: Strong neg. correlations between groups suggests
% the MLE procedure considers these groups similar
% ==================================================
% 
% result = 
%     theta: [0.24235 0.37024 0.3874]
%        SD: [0.040535 0.048536 0.044045]
%       cor: [3x3 double]
%      grpM: [195x1 double]
%        PP: [195x3 double]
%       grp: [1 2 3]
%       txt: {'1'  '2'  '3'}
%     nIter: 4
%       pdf: 'MVN + equal covariance'
%       raw: [0.26154 0.35385 0.38462]


% OUTPUT FROM HISEA.EXE:
% TABLE OF COMPOSITION ESTIMATES
% 
% 
%              RAW      COOK & LORD      COOK        MILLAR      MAXIMUM
%                                   CONSTRAINED  CONSTRAINED  LIKELIHOOD
% 
% 
%  ak1202b   0.2615       0.2570       0.2570       0.2570       0.2365
%  nas1202   0.3538       0.3566       0.3566       0.3566       0.3773
%  ske1202   0.3846       0.3864       0.3864       0.3864       0.3861



% (C) Use Empirical (distribution-free) PDF:
result = f_mle(Y.dat,X.dat,X.stock,[],1,3,100,250)
% 
% ==================================================
%         MAXIMUM LIKELIHOOD ESTIMATION
%             Mixing Proportions:
% --------------------------------------------------
% Group 1 = 0.1497 (0.0381) 
% Group 2 = 0.4119 (0.0509) 
% Group 3 = 0.4384 (0.0445) 
% 
% PDF = Distribution-free (empirical)
% 
% # iterations for EM convergence = 6 
% (Bootstrapped SD's are in parentheses)
% 
% --------------------------------------------------
% Bootstrapped correlations among Groups:
% Group 1 vs. 2 = -0.5310 
% Group 1 vs. 3 = -0.2480 
% Group 2 vs. 3 = -0.6892 
% 
% # iterations for Bootstrap = 250 
% 
% NOTE: Strong neg. correlations between groups suggests
% the MLE procedure considers these groups similar
% ==================================================
% 
% result = 
%     theta: [0.14969 0.41194 0.43837]
%        SD: [0.03806 0.050887 0.044512]
%       cor: [3x3 double]
%      grpM: [195x1 double]
%        PP: [195x3 double]
%       grp: [1 2 3]
%       txt: {'1'  '2'  '3'}
%     nIter: 6
%       pdf: 'Distribution-free (empirical)'
%       raw: [0.16923 0.39487 0.4359]





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


% MLE of mixture proportions based on OTOLITH data (ASSUME EQUAL COV):
resultCod = f_mle(Y.dat,X.dat,X.grp,X.txt,1,1,100,500)
% 
% ==================================================
%         MAXIMUM LIKELIHOOD ESTIMATION
%             Mixing Proportions:
% --------------------------------------------------
% Group 3 = 0.3687 (0.0969) :4T 
% Group 1 = 0.5077 (0.1137) :3Pn4RS 
% Group 5 = 0.0359 (0.0366) :4Vs 
% Group 2 = 0.0877 (0.0560) :3Ps 
% 
% PDF = MVN + equal covariance
% 
% # iterations for EM convergence = 13 
% (Bootstrapped SD's are in parentheses)
% 
% --------------------------------------------------
% Bootstrapped correlations among Groups:
% Group 3 vs. 1 = -0.8240 
% Group 3 vs. 5 = -0.1351 
% Group 3 vs. 2 = +0.0308 
% Group 1 vs. 5 = -0.1691 
% Group 1 vs. 2 = -0.4939 
% Group 5 vs. 2 = -0.0766 
% 
% # iterations for Bootstrap = 500 
% 
% NOTE: Strong neg. correlations between groups suggests
% the MLE procedure considers these groups similar
% ==================================================
% 
% resultCod = 
%     theta: [0.36874 0.50768 0.035916 0.087663]
%        SD: [0.096878 0.11365 0.0366 0.055954]
%       cor: [4x4 double]
%      grpM: [100x1 double]
%        PP: [100x4 double]
%       grp: [3 1 5 2]
%       txt: {'4T'  '3Pn4RS'  '4Vs'  '3Ps'}
%     nIter: 13
%       pdf: 'MVN + equal covariance'
%       raw: [0.33 0.35 0.12 0.2]

% Compare to results with Campana's ISMA.scc Splus code:
% 
% > source("go_ISMA.R")
%         4T     3Pn4RS        4Vs        3Ps 
% 0.37031274 0.50376413 0.03701455 0.08890858 
% [1] 4T     3Pn4RS 4Vs    3Ps   
% Levels: 3Pn4RS 3Ps 4T 4Vs
% [1] 13
% 
% -> actually, results were exactly the same as ISMA.ssc, but their code
%    reports the results of the 2nd to last EM iteration, presumably since
%    it considers the algoritm had already converged.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/iris.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Fisher's Iris data:
load iris.mat

% MLE of mixture proportions (UNEQUAL COV):
resultIris = f_mle(iris,iris,grps,[],1,2,100,250)
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
