% Example Multiple Regression
% 
% by David L. Jones, Jan-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% SEE ALSO: exampleStepwise

% 1) This example uses data from Sokal & Rohlf (1999) Table 16.1. There is 1
% response variable for air pollution (levels of SO2) and 6 explanatory
% variables for 41 US cities.

% Load the file
load SR_table16p1.mat

% Multiple linear regression
model = f_mregress(x,y,1000); 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.66954       0.61122       11.48118      0.00000       0.00100       
% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 111.75344     2.36204       0.02313       0.01400       
%             1 -1.26832      -2.04210      0.04777       0.09600       
%             2 0.06492       4.12300       0.00018       0.00200       
%             3 -0.03928      -2.59613      0.01312       0.01400       
%             4 -3.18162      -1.75304      0.08726       0.08000       
%             5 0.51268       1.41352       0.16524       0.19400       
%             6 -0.05217      -0.32204      0.74910       0.73000       
% ---------------------------------------------------------------------
% 
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================


% 2) This example uses data on temperature (temp) and fish length (len) to
% perform a simple linear regression.

% Load the file:
load USF_regression.mat

% Linear Regression:
model = f_mregress(temp,len,1000); 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.68803       0.68003       86.01357      0.00000       0.00100       
% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 35.30670      2.94759       0.00532       0.00200       
%             1 3.44475       9.27435       0.00000       0.00200       
% ---------------------------------------------------------------------
% 
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================


% 3) This example examines the effect of X ( = [temperature depth current] on
% fish length (Y).

% Load the file:
load USF_multRegress.mat

% Multiple Linear Regression:
model = f_mregress(X,Y,1000); 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.99782       0.99764       5632.68854    0.00000       0.00100       
% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 7.25065       7.38618       0.00000       0.00200       
%             1 2.95450       69.54674      0.00000       0.00200       
%             2 -0.47864      -12.10086     0.00000       0.00200       
%             3 3.97420       86.40984      0.00000       0.00200       
% ---------------------------------------------------------------------
% 
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================


% -----AIC-based forward selection:-----
f_rdaAIC(Y,X,0,1,2,X_txt);
% 
% Performing MARGINAL tests on 3 variables...
% 
% ==================================================
% AIC-based stepwise forward selection (RDA)         
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'           'R2'          'R2adj'       'AIC'       'wts'       'delta'         'ratio'         'var'  
%     [     42765]    [ 0.70999]    [ 0.70255]    [289.26]    [     0]    [   0.99983]    [         1]    'curr' 
%     [     65370]    [  0.5567]    [ 0.54533]    [306.66]    [17.398]    [0.00016674]    [    5996.4]    'temp' 
%     [1.4746e+05]    [     NaN]    [     NaN]    [ 337.8]    [48.538]    [2.8839e-11]    [3.4669e+10]    'none' 
%     [1.4104e+05]    [0.043558]    [0.019034]    [338.19]    [48.926]    [2.3762e-11]    [4.2077e+10]    'depth'
% 
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'         'R2adj'      'AIC'       'wts'        'deltaN'    'var'      'idx'
%     [ 42765]    [0.70999]    [0.70255]    [289.26]    [0.99983]    [48.538]    'curr'     [  3]
%     [1597.2]    [0.98917]    [ 0.9886]    [156.81]    [      1]    [132.45]    'temp'     [  1]
%     [322.18]    [0.99782]    [0.99764]    [93.634]    [      1]    [63.175]    'depth'    [  2]
% 
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

% -> the variables selected are [curr temp depth]


% -----F-ratio based forward selection:-----
f_rdaStepwise(Y,X,[1000],1,[0.05],1,0,X_txt);
% Performing GLOBAL TEST with 1000 permutations...
% 
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% 
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'         'R2adj'  
%     [5632.7]    [0.001]    [0.99782]    [0.99764]
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'         'p '       'R2'          'R2adj'       'Variable'
%     [48.976]    [0.001]    [  0.5567]    [ 0.54533]    'temp'    
%     [1.7761]    [0.187]    [0.043558]    [0.019034]    'depth'   
%     [95.478]    [0.001]    [ 0.70999]    [ 0.70255]    'curr'    
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [   95.478]    [0.001]    [   0.70999]    [      0.70255]    [  0.70255]    'curr'    
%     [   979.45]    [0.001]    [   0.27918]    [       0.2607]    [   0.9886]    'temp'    
%     [   146.43]    [0.001]    [ 0.0086466]    [    -0.016773]    [  0.99764]    'depth'   
% 
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

