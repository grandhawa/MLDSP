% Examples for Multiple Linear Regression
% 
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/Neter_p241.mat'

% Load example data from Neter et al. (1996) page 241:
load Neter_p241.mat
 
% Multiple Linear Regression:
model_1 = f_mregress(x,y,5000,1,1);
% 
% Permuting the data 4999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.91675       0.90750       99.10350      0.00000       0.00020       
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept -68.85707     -1.14729      0.26480       0.23360       
%             1 1.45456       6.86820       0.00000       0.00040       
%             2 9.36550       2.30453       0.03205       0.03560       
% ---------------------------------------------------------------------
% # permutations of residuals =  4999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================
% 
% -> Results are similar to Neter's SYSTAT output and SAS JMP output


% Redundancy Analysis:
result_1 = f_rda(y,x,[0],5000,1);
% 
% Permuting residuals under a full model 4999 times...
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 99.1035    p    =  0.00020 
% R2 = 0.9167   R2adj =  0.90750 
% No. of permutations = 5000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------


% Load data from Sokal & Rohlf (1999) Table 16.1 (air pollution in 41 cities):
load SR_table16p1.mat

% Multiple Linear Regression:
model_2 = f_mregress(x,y,1000,1,1);
% 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.66954       0.61122       11.48118      0.00000       0.00100       
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 111.75344     2.36204       0.02313       0.01200       
%             1 -1.26832      -2.04210      0.04777       0.05000       
%             2 0.06492       4.12300       0.00018       0.00200       
%             3 -0.03928      -2.59613      0.01312       0.00200       
%             4 -3.18162      -1.75304      0.08726       0.07800       
%             5 0.51268       1.41352       0.16524       0.15200       
%             6 -0.05217      -0.32204      0.74910       0.72200       
% ---------------------------------------------------------------------
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================


% Redundancy Analysis:
result_2 = f_rda(y,x,[0],1000,1);
% 
% Permuting residuals under a full model 999 times...
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 11.4812    p    =  0.00100 
% R2 = 0.6695   R2adj =  0.61122 
% No. of permutations = 1000 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------


% -----AIC-based  forward selection:-----
f_rdaAIC(y,x,0,1,2,x_txt);
% 
% Performing MARGINAL tests on 6 variables...
% 
% ==================================================
% AIC-based stepwise forward selection (RDA)         
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'      'R2'           'R2adj'        'AIC'       'delta'     'wts'           'ratio'     'var' 
%     [12876]    [  0.41573]    [  0.40075]    [240.05]    [     0]    [   0.99339]    [     1]    'Manu'
%     [16665]    [  0.24382]    [  0.22443]    [250.62]    [10.574]    [ 0.0050224]    [197.79]    'Pop' 
%     [17895]    [  0.18801]    [  0.16719]    [253.54]    [13.494]    [ 0.0011667]    [851.47]    'Temp'
%     [19028]    [  0.13658]    [  0.11444]    [256.06]    [16.012]    [0.00033126]    [2998.8]    'Days'
%     [22038]    [      NaN]    [      NaN]    [259.87]    [ 19.82]    [4.9357e-05]    [ 20127]    'none'
%     [21840]    [0.0089663]    [-0.016445]    [261.71]    [21.664]    [1.9631e-05]    [ 50604]    'Wind'
%     [21973]    [0.0029467]    [-0.022619]    [261.96]    [21.912]    [1.7339e-05]    [ 57292]    'Rain'
% 
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'         'R2adj'      'AIC'       'wts'        'deltaN'     'var'     'idx'
%     [ 12876]    [0.41573]    [0.40075]    [240.05]    [0.99339]    [  19.82]    'Manu'    [  2]
%     [9116.6]    [0.58632]    [0.56455]    [228.22]    [0.93989]    [ 11.823]    'Pop'     [  3]
%     [9116.6]    [    NaN]    [    NaN]    [228.22]    [0.22822]    [0.73989]    'none'       []
% 
% --------------------------------------------------
% 
% RSS    = residual sum-of-squares 
% R2     = fraction of total variance explained 
% R2adj  = fraction of adjusted total variance explained 
% AIC    = corrected AIC 
% deltaN = delta associated with NO variable addition 
% wts    = AIC weights 
% var    = variable labels 
% idx    = index to selected variables 
% 
% (Note: RSS, R2, and R2adj in Marginal tests are CUMULATIVE) 
% 
% -> the variables selected are [Manu Pop]


% -----F-ratio based forward selection:-----
f_rdaStepwise(y,x,[1000],1,[0.05],1,0,x_txt);

% Performing GLOBAL TEST with 1000 permutations...
% 
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% 
% Now performing stepwise variable selection... 
% 
% Performing CONDITIONAL TESTS with 1000 permutations...
%   -> examining: Temp 
%   -> examining: Manu 
%   -> examining: Pop 
%   -> examining: Wind 
%   -> examining: Rain 
%   -> examining: Days 
% 
% Performing MARGINAL TESTS with 1000 permutations...
%   -> selecting variable 1 of 6 
%   -> selecting variable 2 of 6 
%   -> selecting variable 3 of 6 
% 
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'         'R2adj'  
%     [11.481]    [0.001]    [0.66954]    [0.61122]
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'          'p '       'R2'           'R2adj'        'Variable'
%     [ 9.0301]    [0.003]    [  0.18801]    [  0.16719]    'Temp'    
%     [  27.75]    [0.001]    [  0.41573]    [  0.40075]    'Manu'    
%     [ 12.575]    [0.008]    [  0.24382]    [  0.22443]    'Pop'     
%     [0.35285]    [0.539]    [0.0089663]    [-0.016445]    'Wind'    
%     [0.11526]    [0.738]    [0.0029467]    [-0.022619]    'Rain'    
%     [ 6.1691]    [0.017]    [  0.13658]    [  0.11444]    'Days'    
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [    27.75]    [0.001]    [   0.41573]    [      0.40075]    [  0.40075]    'Manu'    
%     [    15.67]    [0.001]    [   0.17059]    [      0.14933]    [  0.56455]    'Pop'     
% 
% Used 1000 permutation of residuals under a REDUCED model
% (Only variables that significantly contribute to the model are shown)
% --------------------------------------------------
% 
%         R2    = fraction of total variance explained. 
%         R2adj = adjusted R2. 
%    Partial R2 = fraction of variance explained after effects of
%                 variables already in model have been removed. 
% Partial R2adj = adjusted Partial R2
%     Cum R2adj = cumulative fraction of adjusted total variance explained. 
% --------------------------------------------------
% 
% Variable addition halted due to: ALPHA LEVEL 


% -> the variables selected are [Manu Pop]
