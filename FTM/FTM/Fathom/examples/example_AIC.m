% Examples of AIC-based variable selection
% 
% by David L. Jones, May-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% updated Jul-2012


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Fisher's Iris Data:                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Load example data:
load iris.mat

% AIC-based variable selection:
Y    = f_designMatrix(grps,0);
best = f_rdaAIC(Y,iris,0,1,2);
% 
% Performing MARGINAL tests on 4 variables...
% 
% ==================================================
% AIC-based stepwise forward selection (RDA)         
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'       'R2'         'R2adj'      'AIC'        'wts'       'delta'         'ratio'         'var' 
%     [52.931]    [0.47069]    [0.46711]    [-152.16]    [     0]    [   0.70674]    [         1]    '3'   
%     [53.556]    [0.46444]    [0.46082]    [ -150.4]    [1.7592]    [   0.29326]    [      2.41]    '4'   
%     [69.065]    [0.30935]    [0.30469]    [-112.26]    [39.907]    [ 1.526e-09]    [4.6313e+08]    '1'   
%     [79.961]    [0.20039]    [0.19499]    [-90.283]    [61.881]    [ 2.582e-14]    [2.7371e+13]    '2'   
%     [   100]    [    NaN]    [    NaN]    [-58.793]    [93.371]    [3.7491e-21]    [1.8851e+20]    'none'
% 
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'         'R2adj'      'AIC'        'wts'        'deltaN'    'var'     'idx'
%     [52.931]    [0.47069]    [0.46711]    [-152.16]    [0.70674]    [93.371]    '3'       [  3]
%     [44.005]    [0.55995]    [0.55397]    [-177.79]    [0.99753]    [25.623]    '2'       [  2]
%     [40.504]    [0.59496]    [0.58663]    [-188.11]    [0.98852]    [10.321]    '4'       [  4]
%     [40.504]    [    NaN]    [    NaN]    [-188.11]    [0.70815]    [     0]    'none'       []
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Distance-based RDA:                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Repeat Fisher's Iris Data example using db-RDA approach
% 
% Load example data:
load iris.mat

% Create labels:
iris_txt = {'one' 'two' 'three' 'four'};

% AIC-based variable selection:
Y      = f_designMatrix(grps,0);
bestDB = f_rdaDB_AIC(f_dis(Y,'euc'),size(Y,2),iris,0,1,2,0,iris_txt);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        North Sea Ekofisk Oil Field:                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Use an AIC-based stepwise variable selection procedure with a distance-based
% RDA (db-RDA) approach to build an optimal model depicting the
% relationship between a set of environmenal metrics with the benthic
% community in the North Sea. The file 'ekofisk.mat' consists of North Sea
% Ekofisk oil field  data from Gray et al., (1988) and has the following
% variables:
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

% In order to demonstrate how categorical variables can be incorporated
% into this kind of analysis, we will create a new variable specifying the
% geographic zone for each site defined according to their distance from
% the center of drilling activity.
% 
% Create a variable, 'geo', which indicates membership of sites within
% each of 4 geographic zones, defined as follows: 
% 
% zone 1: 0 to 250 m
% zone 2: greater than 250m, but not exceeding 1 km
% zone 3: greater than 1 km, but not exceeding 3.5 km
% zone 4: greater than 3.5 km
% 
% Create grouping variable using logical indices:
geo = zeros(size(env,1),1) + 4; % initialize
geo(env(:,1)<=3.5)  = 3;
geo(env(:,1)<=1)    = 2;
geo(env(:,1)<=0.25) = 1;

% Dummy code the categorical variable:
geoRx = f_dummy(geo);

% Create a new set of environmenal data containg 'geoRx' instead of 'Dist'
env2     = [geoRx env(:,2:end)];
env2_txt = [{'geoRx_01' 'geoRx_02' 'geoRx_03'} env_txt(2:end)];

% Square-root transform the biotic data to down-weight the influence of the
% more abundant species:
biotic_2 = f_normal(biotic,'2');

% Create a Bray-Curtis dissimilarity matrix:
disBC = f_dis(biotic_2,'bc');

% AIC-baseed stepwise variable selection:
best_model = f_rdaDB_AIC(disBC,size(biotic,2),env2,0,1,2,env2_txt);

% Display the variables comprising the optimal model:
env2_txt(best_model.idx)
% ans = 
% 
%     'Sr'    'Ba'

% NOTE: see 'exampleStepwise.m' for a similar example using the
% f_rdaDB_Stepwise function.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Reproduce example for ortho.AIC from 'spacemakeR for R' manual:        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Load example data from 'spacemakeR for R':
load aic.mat

% AIC-based variable selection:
[model,cond] = f_rdaAIC(data.y,data.x,0,1,2);
% 
% Performing MARGINAL tests on 5 variables...
% 
% ==================================================
% AIC-based stepwise forward selection (RDA)         
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'       'R2'           'R2adj'         'AIC'         'wts'       'delta'         'ratio'         'var' 
%     [15.079]    [  0.67575]    [     0.669]    [ -55.681]    [     0]    [         1]    [         1]    '1'   
%     [42.924]    [ 0.076981]    [  0.057751]    [ -3.3738]    [52.307]    [ 4.381e-12]    [2.2826e+11]    '2'   
%     [43.836]    [ 0.057375]    [  0.037737]    [ -2.3229]    [53.358]    [2.5904e-12]    [3.8604e+11]    '5'   
%     [46.504]    [      NaN]    [       NaN]    [ -1.5405]    [54.141]    [1.7518e-12]    [5.7084e+11]    'none'
%     [45.992]    [ 0.011024]    [-0.0095794]    [0.077186]    [55.758]    [7.8019e-13]    [1.2817e+12]    '4'   
%     [ 46.19]    [0.0067522]    [  -0.01394]    [  0.2927]    [55.974]    [7.0049e-13]    [1.4276e+12]    '3'   
% 
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'         'R2adj'      'AIC'        'wts'        'deltaN'     'var'     'idx'
%     [15.079]    [0.67575]    [  0.669]    [-55.681]    [      1]    [ 54.141]    '1'       [  1]
%     [11.499]    [0.75274]    [0.74221]    [-66.967]    [ 0.8646]    [ 11.286]    '2'       [  2]
%     [8.8307]    [0.81011]    [0.79773]    [-77.801]    [0.98872]    [ 10.834]    '5'       [  5]
%     [8.8307]    [    NaN]    [    NaN]    [-77.801]    [0.33206]    [0.51573]    'none'       []
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


% Results from 'spacemakeR':
% 
% $AICc  = -55.68129 -66.96700 -77.80058 -78.31631 -77.65052
% $AICc0 = -1.540526
% $ord   = 1 2 5 4 3
% $R2    = 0.6757547 0.7527357 0.8101103 0.8211346 0.8278868


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Reproduce example for test.W from 'spacemakeR for R' manual:           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Load example data from 'spacemakeR for R':
load aic.mat

% AIC-based variable selection:
best = f_rdaAIC(faudt,MEM,0,0,2);

% Show results:
% best = 
% 
%        RSS: [3x1 double]
%         R2: [3x1 double]
%      R2adj: [3x1 double]
%        AIC: [3x1 double]
%         wt: [3x1 double]
%     deltaN: [3x1 double]
%        var: {' 3'  ' 1'  'none'}
%        idx: [3 1]


% Results from 'spacemakeR:
% $best$AICc  = -91.488003 -94.917598 -95.489469 -95.905881 -96.339367 -96.716195
% $best$ord   =  3 1 2 15 7 24
% $best$AICc0 = -87.47112
