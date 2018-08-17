% Example Stepwise Variable Selection via Forward Addition
% 
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 1) This example uses data from Sokal & Rohlf (1999) Table 16.1. There is 1
% response variable for air pollution (levels of SO2) and 6 explanatory
% variables for 41 US cities.

% Load the file
load SR_table16p1.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AIC-based forward selection:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% -> the variables selected are [Manu Pop ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     F-ratio based forward selection:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_rdaStepwise(y,x,[1000],1,[0.05],1,0,x_txt);
% 
% Performing GLOBAL TEST with 1000 permutations...
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% 
% Now performing stepwise variable selection... 
% 
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'         'R2adj'  
%     [11.481]    [0.001]    [0.66954]    [0.61122]
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'F'          'p '       'R2'           'R2adj'        'Variable'
%     [ 9.0301]    [0.005]    [  0.18801]    [  0.16719]    'Temp'    
%     [  27.75]    [0.001]    [  0.41573]    [  0.40075]    'Manu'    
%     [ 12.575]    [0.004]    [  0.24382]    [  0.22443]    'Pop'     
%     [0.35285]    [ 0.56]    [0.0089663]    [-0.016445]    'Wind'    
%     [0.11526]    [0.721]    [0.0029467]    [-0.022619]    'Rain'    
%     [ 6.1691]    [ 0.02]    [  0.13658]    [  0.11444]    'Days'    
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% Marginal Tests: (sequential variable addition)
% 'Partial F' 'p'     Partial R2'  'Partial R2adj' 'Cum R2adj'  'Variable'
% [27.75]     [0.001] [0.41573]    [0.40075]       [0.40075]    'Manu'    
% [15.67]     [0.001] [0.17059]    [0.14933]       [0.56455]    'Pop'     
% Used 1000 permutation of residuals under a REDUCED model
% (Only variables that significantly contribute to the model are shown)
% --------------------------------------------------
%         R2    = fraction of total variance explained. 
%         R2adj = adjusted R2. 
%    Partial R2 = fraction of variance explained after effects of
%                 variables already in model have been removed. 
% Partial R2adj = adjusted Partial R2
%     Cum R2adj = cumulative fraction of adjusted total variance explained. 
% --------------------------------------------------
% 
% Variable addition halted due to: ALPHA LEVEL 


% Similar results are obtained by S. Dray's 'forward.sel' function in the 
% 'packfor for R' package, which also implements the recommendations of Blanchet
% et al. (in press):
% 
% > result <- forward.sel(y,x,nperm=999,alpha=0.05)
% Testing variable 1
% Testing variable 2
% Testing variable 3
% Procedure stopped (alpha criteria): pvalue for variable 3 is 0.093000
% (superior to  0.050000)
% > result
%   variables order        R2     R2Cum  AdjR2Cum        F  pval
% 1      Manu     2 0.4157267 0.4157267 0.4007453 27.74959 0.001
% 2       Pop     3 0.1705935 0.5863202 0.5645476 15.67046 0.001



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Distance-based Approach:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Use an F-ratio based stepwise variable selection procedure with a distance-based
% RDA (db-RDA) approach to build a model depicting the relationship between
% a set of environmenal metrics with the benthic community in the North
% Sea. The file 'ekofisk.mat' consists of North Sea Ekofisk oil field  data
% from Gray et al., (1988) and has the following variables:
% 
% biotic     = abundance of 174 species of soft-bottom benthic macrofauna from 39 sites
% biotic_txt = cell array of corresponding column labels
% env        = values of 10 environmental variables from the same 39 sites
% env_txt    = cell array of corresponding column labels:
%       'Dist': distance (km) from center of drilling activity
%        'THC': total hydrocarbon concentration
%      'Redox': redox potential (~ organic matter)
%       '%Mud': percent of sediment that was mud
%    Phi_mean': mean sediment particle size
%         'Ba': sediment Barium concentration (nontoxic, tracer for drilling muds)
%         'Sr': sediment Strontium concentration
%         'Cu': sediment Copper concentration
%         'Pb': sediment Lead concentration
%         'Ni': sediment Nickel concentration
% sites_txt  = cell array of site labels
% 
% 
% Gray, J. S., M. Aschan, M. R. Carr, K. R. Clarke, R. H. Green, T. H.
%  Pearson, R. Rosenberg, & R. M. Warwick. 1988. Analysis of community
%  attributes of the benthic macrofauna of Frierfjord/Langesundfjord and in
%  a mesocosm experiment. Mar. Ecol. Prog. Ser. 46: 151-165.

% Load data:
load ekofisk.mat;

% Square-root transform the biotic data to down-weight the influence of the
% more abundant species:
biotic_2 = f_normal(biotic,'2');

% Create a Bray-Curtis dissimilarity matrix:
disBC = f_dis(biotic_2,'bc');

% Forward selection of variables:
f_rdaDB_Stepwise(disBC,size(biotic,2),env,[1000],1,0.05,1,[0],env_txt);

% Performing GLOBAL TEST with 1000 permutations...
% 
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% 
% Now performing stepwise variable selection... 
% 
% Performing CONDITIONAL TESTS with 1000 permutations...
%   -> examining: Dist 
%   -> examining: THC 
%   -> examining: Redox 
%   -> examining: %Mud 
%   -> examining: Phi_mean 
%   -> examining: Ba 
%   -> examining: Sr 
%   -> examining: Cu 
%   -> examining: Pb 
%   -> examining: Ni 
% 
% Performing MARGINAL TESTS with 1000 permutations...
%   -> selecting variable 1 of 10 
%   -> selecting variable 2 of 10 
% 
% ==================================================
%                 Stepwise db-RDA:                 
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'         'R2adj'  
%     [3.1872]    [0.001]    [0.53233]    [0.36531]
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'          'p '       'R2'          'R2adj'        'Variable'
%     [ 9.6785]    [0.001]    [ 0.20734]    [  0.18592]    'Dist'    
%     [ 5.3816]    [0.001]    [ 0.12698]    [  0.10339]    'THC'     
%     [ 1.1218]    [0.316]    [0.029427]    [0.0031949]    'Redox'   
%     [ 7.3029]    [0.001]    [ 0.16484]    [  0.14227]    '%Mud'    
%     [ 4.8246]    [0.001]    [ 0.11535]    [ 0.091443]    'Phi_mean'
%     [ 6.4462]    [0.001]    [ 0.14837]    [  0.12535]    'Ba'      
%     [ 10.775]    [0.001]    [ 0.22554]    [   0.2046]    'Sr'      
%     [ 5.2993]    [0.001]    [ 0.12528]    [  0.10164]    'Cu'      
%     [ 9.4252]    [0.001]    [ 0.20302]    [  0.18148]    'Pb'      
%     [0.51427]    [0.945]    [0.013709]    [-0.012948]    'Ni'      
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [   10.775]    [0.001]    [   0.22554]    [       0.2046]    [   0.2046]    'Sr'      
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
% Variable addition halted due to: ALPHA LEVEL 


