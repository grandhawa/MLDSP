% Example of a PCNM-based analysis of spatial ecology
% 
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/oribatid_mites.mat'
% 
% Nov-2011: edited to work with updated functions

% -----Notes:-----
% This example follows the analysis of Oribatid mites presented in Borcard
% et al. (2004).

% -----Variables:-----
% bio = structure of abundances of 35 species of mites from 70 cores:
%  .num = # per 5 cm diameter, 7 cm deep core
%  .txt = species codes
% 
% env = structure of environmental variables with the following fields:
%  .x       = x coordinates (meters)
%  .y       = y coordinates (meters)
%  .density = bulk density  (grams/liter of dry matter)
%  .water   = water content (grams/liter of fresh matter)
%  .substra = substrate (1 = Sphagn_1, 2 = Sphagn_2, 3 = Sphagn_3, 4 = Sphagn_4,
%             5 = Lignlitt, 6 = Barepeat, 7 = Interface)   
%  .schrub  = schrub coverage (density classes ranging from 1-3)
%  .micro   = microtopography (1 = blanket, 2 = hummock)

% Load data:
load oribatid_mites.mat

% Hellinger transform the species data:
bio.H = f_hellinger(bio.num);

% Dummy code qualitative variables:
env.Qsubstra = f_dummy(env.substra,1);
env.Qschrub  = f_dummy(env.schrub,1);
env.Qmicro   = f_dummy(env.micro,1);

% Perform preliminary RDA on X,Y coordinates to identify linear trends:
result_XY = f_rda(bio.H,[env.x env.y],[0],1000,1);
% 
% Permuting residuals under a full model 999 times...
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 13.2798    p    =  0.00100 
% R2 = 0.2839   R2adj =  0.26250 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.1058  0.0062
% Residual Eigenvalues:
%   0.0936  0.0354  0.0299  0.0167...
% 
% Species-Environment Correlations (r):
%   0.8308  0.5157
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.2839): 
%   0.2683  0.0156
% Cumulative: 
%   0.2683  0.2839
% 
% Residual axes  (total = 0.7161):
%   0.2373  0.0898  0.0758  0.0423...
% Cumulative: 
%   0.2373  0.3271  0.4030  0.4452...
% =================================================
% 
% -> the X,Y coordinates significantly explained 28.4% (or 26.3% if you use
% R^2adj) of the variance in the species data.


% Perform multiple regression of the 1st RDA axis on the set of environmental
% variables: 
model_XY = f_mregress([env.density env.water env.Qsubstra env.Qschrub env.Qmicro],...
   result_XY.siteScores(:,1),1000,1,1);
% 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.83816       0.80747       27.30783      0.00000       0.00100       
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept -0.51905      -5.47439      0.00000       0.00200       
%             1 -0.01175      -5.22956      0.00000       0.00200       
%             2 0.00174       8.04996       0.00000       0.00200       
%             3 -0.02979      -0.56286      0.57536       0.62600       
%             4 0.06999       1.10208       0.27426       0.30200       
%             5 -0.03642      -0.19856      0.84319       0.89000       
%             6 -0.10974      -0.83494      0.40663       0.48000       
%             7 0.13975       1.05638       0.29448       0.27600       
%             8 0.09868       0.72835       0.46886       0.36800       
%             9 0.15754       2.18200       0.03252       0.03400       
%            10 0.13043       2.52307       0.01394       0.00200       
%            11 0.27252       4.77407       0.00001       0.00200       
% ---------------------------------------------------------------------
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================
% 
% -> Note that the 1st canonical axis can be significantly explained by density,
% water, schrub, and microhabitat (but not substratum).

% Create PCNM's:
pcnm     = f_pcnm([env.x env.y]);
pcnm.txt = cellstr(num2str([1:size(pcnm.evects,2)]')); % create text labels

% Create PCNM's (alternate method):
mst = f_mst(f_dis([env.x env.y],'euc'),[env.x env.y]); % Minimum Spanning Tree
A   = f_dis2sim(mst.tDis,2,2);                         % weighting matrix
MEM = f_eigenMaps(mst,A,1000,0);                       % Moran's eigenvector maps

% Use PCNM's to model the spatial structure of the biotic data; since the
% preliminary RDA indicated a significant linear gradient in the data, use
% linearly de-trended data for subsequent analyses (i.e., the residuals from the
% preliminary analysis to perform a partial RDA):
%
% Residuals serve as the de-trended data:
bio.dt = result_XY.res;
 
% Select only significant PCNM's:
f_rdaStepwise(bio.dt,pcnm.evects,[1000 0],1,[0.05]);
% 
% Performing GLOBAL TEST with 1000 permutations...
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% Now performing stepwise variable selection... 
% Performing CONDITIONAL TESTS with 0 permutations...
% Performing MARGINAL TESTS with 1000 permutations...
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'        'R2adj'  
%     [1.6549]    [0.001]    [0.7324]    [0.28983]
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'          'p '       'R2'           'R2adj'          'Variable'
%     [ 4.7304]    [0.002]    [ 0.065041]    [   0.051291]    ' 1'      
%     [0.98282]    [  NaN]    [ 0.014247]    [-0.00024908]    ' 2'      
%     [ 3.2983]    [  NaN]    [  0.04626]    [   0.032235]    ' 3'      
%     [ 3.7439]    [  NaN]    [ 0.052184]    [   0.038246]    ' 4'      
%     [ 1.6374]    [  NaN]    [ 0.023513]    [  0.0091531]    ' 5'      
%     [ 2.5954]    [  NaN]    [ 0.036765]    [     0.0226]    ' 6'      
%     [ 2.6245]    [  NaN]    [ 0.037161]    [   0.023001]    ' 7'      
%     [0.78298]    [  NaN]    [ 0.011383]    [ -0.0031552]    ' 8'      
%     [0.78959]    [  NaN]    [ 0.011478]    [ -0.0030588]    ' 9'      
%     [ 1.9232]    [  NaN]    [ 0.027504]    [   0.013203]    '10'      
%     [ 2.9229]    [  NaN]    [ 0.041213]    [   0.027113]    '11'      
%     [0.45607]    [  NaN]    [0.0066622]    [ -0.0079458]    '12'      
%     [  1.058]    [  NaN]    [ 0.015321]    [ 0.00084032]    '13'      
%     [0.75302]    [  NaN]    [ 0.010952]    [ -0.0035923]    '14'      
%     [0.51337]    [  NaN]    [ 0.007493]    [ -0.0071027]    '15'      
%     [0.98883]    [  NaN]    [ 0.014333]    [-0.00016197]    '16'      
%     [0.34694]    [  NaN]    [0.0050761]    [ -0.0095551]    '17'      
%     [0.71019]    [  NaN]    [ 0.010336]    [ -0.0042179]    '18'      
%     [0.69444]    [  NaN]    [ 0.010109]    [ -0.0044481]    '19'      
%     [ 1.6426]    [  NaN]    [ 0.023586]    [  0.0092265]    '20'      
%     [0.46416]    [  NaN]    [0.0067796]    [ -0.0078265]    '21'      
%     [0.59183]    [  NaN]    [0.0086283]    [ -0.0059507]    '22'      
%     [ 1.6491]    [  NaN]    [ 0.023678]    [  0.0093202]    '23'      
%     [0.88022]    [  NaN]    [ 0.012779]    [  -0.001739]    '24'      
%     [0.70847]    [  NaN]    [ 0.010311]    [  -0.004243]    '25'      
%     [0.57493]    [  NaN]    [ 0.008384]    [ -0.0061986]    '26'      
%     [0.62071]    [  NaN]    [0.0090455]    [ -0.0055274]    '27'      
%     [0.78539]    [  NaN]    [ 0.011418]    [ -0.0031199]    '28'      
%     [0.96601]    [  NaN]    [ 0.014007]    [-0.00049279]    '29'      
%     [0.71107]    [  NaN]    [ 0.010349]    [ -0.0042049]    '30'      
%     [0.46503]    [  NaN]    [0.0067923]    [ -0.0078137]    '31'      
%     [0.61044]    [  NaN]    [0.0088972]    [ -0.0056778]    '32'      
%     [0.92197]    [  NaN]    [ 0.013377]    [ -0.0011321]    '33'      
%     [0.64958]    [  NaN]    [0.0094622]    [ -0.0051045]    '34'      
%     [0.72647]    [  NaN]    [  0.01057]    [   -0.00398]    '35'      
%     [0.94124]    [  NaN]    [ 0.013653]    [-0.00085232]    '36'      
%     [ 2.5494]    [  NaN]    [ 0.036137]    [   0.021962]    '37'      
%     [0.55218]    [  NaN]    [0.0080549]    [ -0.0065325]    '38'      
%     [ 1.0275]    [  NaN]    [ 0.014885]    [ 0.00039842]    '39'      
%     [0.61884]    [  NaN]    [0.0090185]    [ -0.0055547]    '40'      
%     [0.32101]    [  NaN]    [0.0046986]    [ -0.0099382]    '41'      
%     [0.40642]    [  NaN]    [0.0059412]    [ -0.0086773]    '42'      
%     [0.33587]    [  NaN]    [ 0.004915]    [ -0.0097186]    '43'      
% 
% 
% Used 0 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [   4.7304]    [0.002]    [  0.065041]    [     0.051291]    [ 0.051291]    ' 1'      
%     [   3.9606]    [0.002]    [  0.052184]    [     0.038246]    [ 0.090873]    ' 4'      
%     [   3.6499]    [0.002]    [   0.04626]    [     0.032235]    [  0.12546]    ' 3'      
%     [   3.3683]    [0.006]    [  0.041213]    [     0.027113]    [  0.15576]    '11'      
%     [    3.137]    [0.005]    [  0.037161]    [     0.023001]    [  0.18263]    ' 7'      
%     [   3.2108]    [0.002]    [  0.036765]    [       0.0226]    [  0.20992]    ' 6'      
%     [   3.2696]    [0.005]    [  0.036137]    [     0.021962]    [  0.23739]    '37'      
%     [   2.5508]    [0.012]    [  0.027504]    [     0.013203]    [    0.256]    '10'      
%     [   2.2406]    [0.014]    [  0.023678]    [    0.0093202]    [  0.27083]    '23'      
%     [   2.2795]    [ 0.02]    [  0.023586]    [    0.0092265]    [  0.28606]    '20'      
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
% 
% Variable addition halted due to: R^2adj 
% 
% -> Above are the top 10 PCNM's. A total of 12 PCNM were selected in the
% analysis of Borcard et al. (2004); however, they had more liberal stopping
% rules. They report that their subset of PCNM's explained 45.1% of the
% detrended data, but this would have been based on R^2 and not R^2adj, the
% latter being a more accurate measure.

% Create index to selected PCNM's:
pcnm.sel = [1 4 3 11 7 6 37 10 23 20]; % index to selected PCNM's

% Re-do RDA using the de-trended data and the selected PCNM's:
result_PCNM = f_rda(bio.dt,pcnm.evects(:,pcnm.sel),[0],1000,1);
% 
% Permuting residuals under a full model 999 times...
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 3.7647    p    =  0.00100 
% R2 = 0.3895   R2adj =  0.28606 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.0595  0.0220  0.0099  0.0067  0.0049  0.0024  0.0017  0.0013  0.0009  0.0006
% Residual Eigenvalues:
%   0.0372  0.0235  0.0180  0.0121  0.0104  0.0098  0.0077  0.0063  0.0060 0.0047...
% 
% Species-Environment Correlations (r):
%   0.8042  0.8174  0.6149  0.6985  0.6253  0.5366  0.4981  0.4024  0.5092  0.3206
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.3895): 
%   0.2107  0.0778  0.0351  0.0239  0.0174  0.0084  0.0061  0.0046  0.0033  0.0023
% Cumulative: 
%   0.2107  0.2886  0.3236  0.3475  0.3649  0.3733  0.3794  0.3840  0.3872  0.3895
% 
% Residual axes  (total = 0.6105):
%   0.1317  0.0831  0.0636  0.0428  0.0368  0.0348  0.0273  0.0224  0.0211 0.0168 ...
% Cumulative: 
%   0.1317  0.2147  0.2784  0.3212  0.3580  0.3927  0.4200  0.4424  0.4635 0.4803 ...
% ==================================================

% Create canonical plot to see which PCNM's are most important for each axis:
f_rdaPlot(result_PCNM,bio.H,1,[3],[0.02],bio.txt,pcnm.txt(pcnm.sel),'none');
% 
% Examine the PCNM vectors and note that Axis I, which explains 21.07% of the
% de-trended data, is mainly a combination of PCNM 4, 11, 7, 3, and 1 (in
% decreasing order of importance). Axis II, which explains 7.78% of the
% de-trended data, is mainly a combination of PCNM 1, 6, 3, and 20. 

%-----Interpret Canonical Axes I & II:-----
% 
% Aggregate environmental variables:
X.dat = [env.density env.water env.Qsubstra env.Qschrub env.Qmicro];
% 
% 1) Use Multiple Regression approach (Borcard et al., 2004) to interpret
% canonical Axis I in terms of the environmental variables: 
% 
f_mregress(X.dat,result_PCNM.siteScores(:,1),1000,1,1);
% 
% Permuting the data 999 times...
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.48122       0.38284       4.89105       0.00002       0.00100       
% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 0.58534       4.45044       0.00003       0.00200       
%             1 -0.00449      -1.44066      0.15420       0.25000       
%             2 -0.00056      -1.87857      0.06453       0.11000       
%             3 -0.01092      -0.14871      0.88222       0.78200       
%             4 -0.06486      -0.73619      0.46411       0.21600       
%             5 -0.20075      -0.78891      0.43286       0.53200       
%             6 0.22744       1.24749       0.21644       0.25000       
%             7 -0.47574      -2.59238      0.01163       0.00200       
%             8 -0.09108      -0.48462      0.62948       0.78400       
%             9 -0.09460      -0.94458      0.34817       0.32200       
%            10 -0.22788      -3.17772      0.00222       0.00200       
%            11 -0.06472      -0.81740      0.41651       0.29000       
% ---------------------------------------------------------------------
% 
% # permutations of residuals =   999 
% F-test is one-tailed, t-tests are two-tailed 
% =====================================================================
% 
% 
% -> variable 10 is the most significant, and corresponds to one of the Qschrub
% variables; note the multiple R^2 = 0.485.
% 
% 
% 2) Use superior Stepwise RDA approach to interpret canonical Axis I in terms
% of the environmental variables:  
% 
% Group categorical variables, aggregate their labels:
X.grp = [1 2 3 3 3 3 3 3 4 4 5]; 
X.txt = {'density' 'water' 'Qsubstra' 'Qschrub' 'Qmicro'};
% 
f_rdaStepwise(result_PCNM.siteScores(:,1),X.dat,1000,1,0.05,1,X.grp,X.txt);
% 
% Performing GLOBAL TEST with 1000 permutations...
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% Now performing stepwise variable selection... 
% Performing CONDITIONAL TESTS with 1000 permutations...
% Performing MARGINAL TESTS with 1000 permutations...
% 
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'        'p '       'R2'         'R2adj'  
%     [4.891]    [0.001]    [0.48122]    [0.38284]
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'         'p '       'R2'          'R2adj'       'Variable'
%     [ 14.68]    [0.001]    [ 0.17755]    [ 0.16545]    'density' 
%     [18.005]    [0.001]    [ 0.20935]    [ 0.19772]    'water'   
%     [2.4884]    [0.031]    [ 0.19158]    [ 0.11459]    'Qsubstra'
%     [ 8.688]    [0.001]    [ 0.20594]    [ 0.18223]    'Qschrub' 
%     [6.8474]    [0.012]    [0.091485]    [0.078125]    'Qmicro'  
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [   18.005]    [0.001]    [   0.20935]    [      0.19772]    [  0.19772]    'water'   
%     [   7.2317]    [0.011]    [  0.077026]    [     0.063453]    [  0.26507]    'density' 
%     [   4.9773]    [0.014]    [  0.094775]    [     0.067754]    [  0.34307]    'Qschrub' 
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
% 
% -> More information is provided here than in the multiple regression shown
% above. Canonical Axis I can be interpreted in terms of water, density, and
% schrubs, which account for a combined  R^2 = 0.387 (R^2adj = 0.349).
% 
% 
% Do the same for AXIS II:
f_rdaStepwise(result_PCNM.siteScores(:,2),X.dat,1000,1,0.05,1,X.grp,X.txt);
% 
% Performing GLOBAL TEST with 1000 permutations...
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% Now performing stepwise variable selection... 
% Performing CONDITIONAL TESTS with 1000 permutations...
% Performing MARGINAL TESTS with 1000 permutations...
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'         'R2adj'  
%     [4.5542]    [0.001]    [0.46344]    [0.36168]
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'          'p '       'R2'          'R2adj'          'Variable'
%     [0.98265]    [0.321]    [0.014245]    [-0.00025155]    'density' 
%     [ 1.6653]    [0.203]    [0.023904]    [  0.0095498]    'water'   
%     [ 1.1887]    [0.339]    [  0.1017]    [   0.016146]    'Qsubstra'
%     [  4.629]    [0.009]    [  0.1214]    [   0.095178]    'Qschrub' 
%     [ 4.0098]    [0.057]    [0.055684]    [   0.041797]    'Qmicro'  
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [    4.629]    [0.009]    [    0.1214]    [     0.095178]    [ 0.095178]    'Qschrub' 
%     [   9.2749]    [0.004]    [   0.10825]    [     0.095141]    [  0.19464]    'Qmicro'  
%     [   8.8344]    [0.001]    [  0.092172]    [     0.078822]    [   0.2801]    'density' 
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
% 
% -> Canonical Axis II can be interpreted in terms of schrubs, microtopography, and
% density, which account for a combined R^2adj = 0.28. Borcard et al.'s (2004)
% multiple regression technique was unable to say much about this axis and
% forced them to conclude that the patterns were produced by some unmeasured
% factors.
