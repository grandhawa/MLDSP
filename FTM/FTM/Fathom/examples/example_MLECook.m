% Example for maximum likelihood estimation via Cook's constrained
% corrected classification
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
result = f_mleCook(Y.dat,X.dat,X.stock,[],1,1);
% ========================================================
%  COOK'S CONSTRAINED CORRECTED CLASSIFICATION ESTIMATOR 
%             Mixing Proportions:
% --------------------------------------------------------
% Group 1 = 0.2576 (NaN) 
% Group 2 = 0.3507 (NaN) 
% Group 3 = 0.3917 (NaN) 
% 
% PDF = MVN + equal covariance
% (Bootstrapped SD's are in parentheses)
% 
% ==================================================

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
resultCod = f_mleCook(Y.dat,X.dat,X.grp,X.txt,1,1);
% 
% ========================================================
%  COOK'S CONSTRAINED CORRECTED CLASSIFICATION ESTIMATOR 
%             Mixing Proportions:
% --------------------------------------------------------
% Group 3 = 0.2441 (NaN) :4T 
% Group 1 = 0.6145 (NaN) :3Pn4RS 
% Group 5 = 0.0492 (NaN) :4Vs 
% Group 2 = 0.0922 (NaN) :3Ps 
% 
% PDF = MVN + equal covariance
% (Bootstrapped SD's are in parentheses)
% 
% ==================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/iris.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Fisher's Iris data:
load iris.mat

% MLE of mixture proportions (UNEQUAL COV):
resultIris = f_mleCook(iris,iris,grps,[],1,2);
% 
% ========================================================
%  COOK'S CONSTRAINED CORRECTED CLASSIFICATION ESTIMATOR 
%             Mixing Proportions:
% --------------------------------------------------------
% Group 1 = 0.3333 (NaN) 
% Group 2 = 0.3406 (NaN) 
% Group 3 = 0.3261 (NaN) 
% 
% PDF = MVN + unequal covariances
% (Bootstrapped SD's are in parentheses)
% 
% ==================================================








