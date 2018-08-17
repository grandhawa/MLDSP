% Example of using RDA with covariables
% 
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/covar.mat'

% -----References:-----
% Follow along in pages 528-531 in Legendre & Legendre (1998). The data is from
% Table 10.5.

% Load data file:
load covar.mat

% Regress Y against XW using Multiple Linear Regression:
model = f_mregress([x w],y,1000,1,1);
% 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.58354       0.39133       3.03596       0.03997       0.03900       
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 9.63671       10.45185      0.00000       0.00200       
%             1 -0.89572      -2.04450      0.05501       0.05400       
%             2 -1.33929      -1.64987      0.11541       0.09200       
%             3 0.53663       0.47113       0.64291       0.52600       
%             4 0.10247       1.24795       0.22721       0.19400       
%             5 0.13503       1.19856       0.24543       0.25200       
%             6 0.01784       1.43550       0.16740       0.11000       
% ---------------------------------------------------------------------
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================


% Regress Y on XW using RDA:
result = f_rda(y,[x w],[0],1000,1);
% 
% Permuting residuals under a full model 999 times...
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 3.0360    p    =  0.06200 
% R2 = 0.5835   R2adj =  0.39133 
% No. of permutations = 1000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------




% -----Y vs. X after the effect of W removed:-----

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            LONG METHOD:                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) Remove effect of W on both Y & X:
result_Y = f_rda(y,w,[0],1000,0);
result_X = f_rda(x,w,[0],1000,0);
% 
% 2) Regression of Yres on Xres:
result_Res = f_rda(result_Y.res,result_X.res,[0],1000,1);
% 
% Permuting residuals under a full model 999 times...
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 2.5062    p    =  0.08300 
% R2 = 0.3197   R2adj =  0.19213 
% No. of permutations = 1000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------
% 
% 3) Regression of Y on Xres:
result_Res2 = f_rda(y,result_X.res,[0],1000,1);
% 
% Permuting residuals under a full model 999 times...
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 1.2977    p    =  0.23100 
% R2 = 0.1957   R2adj =  0.04490 
% No. of permutations = 1000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------
% -> this R2 directly provides the proportion of the variation in Y accounted
% for by X after the effects of W are removed.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            SHORT METHOD:                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
result_cov = f_rda(y,x,w,1000,1);
% 
% Permuting residuals under a reduced model 999 times...
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 2.0363    p    =  0.15400 
% R2 = 0.1957   R2adj =  0.04490 
% No. of permutations = 1000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------
% -> this produces the same R2 as step 3 above, but the F-ratios are different


% Try with db-RDA:
yDis          = f_dis(y,'euc');
result_cov_DB = f_rdaDB(yDis,size(y,2),x,w,1000,1);

% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 2.0363    p    =  0.27800 
% R2 = 0.1957   R2adj =  0.04490 
% No. of permutations = 1000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------
