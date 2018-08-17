% Example Random Forest Regression 2
% 
% -----Author:-----
% by David L. Jones, Oct-2011
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% File: '.../examples/boston_housing.csv'


% Import Boston Housing data:
raw = f_importCSV('boston_housing.csv',1);
% 
% Variables:
% medv    = median value of owner-occupied homes in $1000's
% crim    = per capita crime rate by town
% zn      = proportion of residential land zoned for lots over 25,000 sq.ft.
% indus   = proportion of non-retail business acres per town
% chas    = Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)
% nox     = nitric oxides concentration (parts per 10 million)
% rm      = average number of rooms per dwelling
% age     = proportion of owner-occupied units built prior to 1940
% dis     = weighted distances to five Boston employment centres
% rad     = index of accessibility to radial highways
% tax     = full-value property-tax rate per $10,000
% ptratio = pupil-teacher ratio by town
% b       = 1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town
% lstat   = % lower status of the population
% lon     = longitude
% lat     = latitude
% 
% These are Boston house-price data from: Harrison, D. and Rubinfeld, D. L.
% 1978. Hedonic prices and the demand for clean air. J. Environ. Economics
% & Management 5: 81-102.

% Parse the data:
Y     = raw.dat(:,1);                           % response variable (medv)
[e,n] = f_deg2utm(raw.dat(:,16),raw.dat(:,15)); % convert lat/lon to UTM coordinates
X     = [raw.dat(:,2:14) e n];                  % predictors
X_txt = [raw.txt(2:14) {'e' 'n'}];              % variable labels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            RANDOM FOREST:              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a RF regression:
rf = f_RFreg(X,Y,[],[],1,0,'stnd',1);

% Show statistics:
[rf.mse(end) rf.rsq(end)]
% ans = 9.3191      0.88961

% Create diagnostic plots:
f_RFregPlot(rf,1,1,1);
% -> variables 6 & 13 seem to be the most important
% X_txt([6 13]) = 'rm' 'lstat'

% Perform AIC-based variable selection:
f_RFregAIC(X,Y,1000,0,1,2,X_txt);

% ==================================================
% AIC-based stepwise forward selection (RANDOM FOREST)
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'      'R2'           'AIC'       'wts'       'delta'          'ratio'          'var'    
%     [18494]    [  0.56704]    [  1825]    [     0]    [    0.99799]    [          1]    'lstat'  
%     [18954]    [  0.55629]    [1837.4]    [12.419]    [  0.0020066]    [     497.35]    'nox'    
%     [23004]    [  0.46148]    [1935.4]    [110.41]    [ 1.0582e-24]    [ 9.4308e+23]    'rm'     
%     [24315]    [  0.43077]    [1963.4]    [138.47]    [ 8.5382e-31]    [ 1.1689e+30]    'indus'  
%     [28602]    [  0.33042]    [2045.6]    [220.62]    [ 1.2347e-48]    [ 8.0827e+47]    'ptratio'
%     [28988]    [  0.32137]    [2052.4]    [227.41]    [ 4.1374e-50]    [ 2.4121e+49]    'tax'    
%     [34054]    [  0.20278]    [2133.9]    [308.91]    [ 8.3145e-68]    [ 1.2003e+67]    'rad'    
%     [34412]    [   0.1944]    [2139.2]    [ 314.2]    [ 5.9012e-69]    [ 1.6912e+68]    'zn'     
%     [40379]    [ 0.054713]    [2220.1]    [395.11]    [  1.591e-86]    [ 6.2727e+85]    'crim'   
%     [41475]    [ 0.029067]    [2233.6]    [408.66]    [ 1.8213e-89]    [ 5.4797e+88]    'e'      
%     [41852]    [  0.02023]    [2238.2]    [413.24]    [   1.84e-90]    [ 5.4239e+89]    'chas'   
%     [42888]    [      NaN]    [2248.6]    [ 423.6]    [ 1.0372e-92]    [ 9.6221e+91]    'none'   
%     [44129]    [-0.033075]    [  2265]    [440.05]    [ 2.7787e-96]    [ 3.5916e+95]    'dis'    
%     [45040]    [-0.054397]    [2275.3]    [450.39]    [ 1.5819e-98]    [ 6.3089e+97]    'b'      
%     [50937]    [ -0.19245]    [2337.6]    [512.64]    [4.7879e-112]    [2.0844e+111]    'age'    
%     [58135]    [ -0.36095]    [2404.5]    [579.52]    [1.4368e-126]    [6.9457e+125]    'n'      
% 
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'         'AIC'       'wts'        'deltaN'    'var'      'idx'
%     [ 18494]    [0.56704]    [  1825]    [0.99799]    [ 423.6]    'lstat'    [ 13]
%     [9910.1]    [  0.768]    [1511.3]    [0.99428]    [313.67]    'indus'    [  3]
%     [7552.2]    [ 0.8232]    [1375.8]    [      1]    [135.46]    'rm'       [  6]
%     [6150.2]    [0.85602]    [  1274]    [      1]    [101.87]    'nox'      [  5]
%     [6093.8]    [0.85734]    [1271.3]    [0.44184]    [  2.61]    'dis'      [  8]
%     [  4983]    [0.88335]    [1171.6]    [0.92878]    [99.772]    'e'        [ 14]
%     [  4983]    [    NaN]    [1171.6]    [0.31393]    [1.5173]    'none'        []
% 
% --------------------------------------------------
% 
% # Trees = 1000
% RSS     = residual sum of squares 
% R2      = fraction of total variation explained 
% AIC     = corrected AIC 
% deltaN  = delta associated with NO variable addition 
% wts     = AIC weights 
% var     = variable labels 
% idx     = index to selected variables 
% 
% (Note: RSS, R2 in Marginal Tests are CUMULATIVE) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             RDA + AIC               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
f_rdaAIC(Y,X,0,1,2,X_txt);
% 
% ==================================================
% AIC-based stepwise forward selection (RDA)         
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'      'R2'            'R2adj'        'AIC'       'wts'       'delta'         'ratio'         'var'    
%     [19472]    [   0.54415]    [  0.54324]    [  1851]    [     0]    [         1]    [         1]    'lstat'  
%     [22062]    [   0.48353]    [   0.4825]    [1914.2]    [63.176]    [1.9119e-14]    [5.2303e+13]    'rm'     
%     [31702]    [   0.25785]    [  0.25637]    [2097.6]    [246.62]    [2.8058e-54]    [3.5641e+53]    'ptratio'
%     [32721]    [   0.23399]    [  0.23247]    [2113.7]    [262.63]    [9.3657e-58]    [1.0677e+57]    'indus'  
%     [33339]    [   0.21953]    [  0.21798]    [2123.1]    [272.09]    [8.2443e-60]    [ 1.213e+59]    'tax'    
%     [34916]    [    0.1826]    [  0.18098]    [2146.5]    [295.48]    [6.8753e-65]    [1.4545e+64]    'nox'    
%     [36276]    [   0.15078]    [   0.1491]    [2165.8]    [314.81]    [4.3732e-69]    [2.2867e+68]    'crim'   
%     [36495]    [   0.14564]    [  0.14394]    [2168.9]    [317.86]    [9.4955e-70]    [1.0531e+69]    'rad'    
%     [36647]    [   0.14209]    [  0.14039]    [  2171]    [319.95]    [3.3319e-70]    [3.0012e+69]    'age'    
%     [37167]    [   0.12992]    [  0.12819]    [2178.1]    [327.08]    [9.4293e-72]    [1.0605e+71]    'zn'     
%     [37966]    [    0.1112]    [  0.10943]    [2188.9]    [337.86]    [4.3146e-74]    [2.3177e+73]    'b'      
%     [38369]    [   0.10178]    [ 0.099998]    [2194.2]    [343.19]    [2.9995e-75]    [3.3339e+74]    'e'      
%     [40048]    [  0.062464]    [ 0.060604]    [2215.9]    [364.87]    [5.8876e-80]    [1.6985e+79]    'dis'    
%     [41404]    [  0.030716]    [ 0.028793]    [2232.8]    [381.72]    [1.2905e-83]    [7.7491e+82]    'chas'   
%     [42716]    [       NaN]    [      NaN]    [2246.5]    [395.49]    [1.3201e-86]    [7.5753e+85]    'none'   
%     [42704]    [0.00028455]    [-0.001699]    [2248.4]    [397.36]    [5.1775e-87]    [1.9314e+86]    'n'      
% 
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'      'R2'         'R2adj'      'AIC'       'wts'        'deltaN'    'var'        'idx'
%     [19472]    [0.54415]    [0.54324]    [  1851]    [      1]    [395.49]    'lstat'      [ 13]
%     [15439]    [0.63856]    [0.63712]    [1735.6]    [      1]    [115.41]    'rm'         [  6]
%     [13728]    [0.67862]    [ 0.6767]    [1678.2]    [      1]    [57.413]    'ptratio'    [ 11]
%     [13229]    [0.69031]    [0.68784]    [1661.5]    [0.80593]    [16.698]    'dis'        [  8]
%     [12469]    [0.70809]    [0.70517]    [1633.6]    [0.99411]    [27.872]    'nox'        [  5]
%     [12141]    [0.71577]    [0.71236]    [1622.2]    [ 0.5632]    [11.443]    'chas'       [  4]
%     [11868]    [0.72216]    [0.71826]    [1612.8]    [0.84716]    [9.4358]    'b'          [ 12]
%     [11678]    [0.72661]    [0.72221]    [1606.7]    [0.61512]    [6.0902]    'zn'         [  2]
%     [11584]    [0.72883]    [ 0.7239]    [1604.6]    [0.31192]    [2.0389]    'crim'       [  1]
%     [11355]    [0.73418]    [0.72881]    [1596.6]    [0.94127]    [7.9959]    'rad'        [  9]
%     [11081]    [0.74058]    [0.73481]    [1586.4]    [0.98093]    [10.244]    'tax'        [ 10]
%     [11081]    [    NaN]    [    NaN]    [1586.4]    [ 0.3313]    [     0]    'none'          []
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     F-ratio based forward selection:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
f_rdaStepwise(Y,X,1000,1,0.05,1,0,X_txt);
% 
% Performing GLOBAL TEST with 1000 permutations...
% 
% GLOBAL TEST (p = 0.001) is significant at alpha = 0.05
% 
% Now performing stepwise variable selection... 
% 
% ==================================================
% Stepwise REDUNDANCY ANALYSIS:
% --------------------------------------------------
% Global Test: (all variables included)
%     'F'         'p '       'R2'         'R2adj'  
%     [93.976]    [0.001]    [0.74206]    [0.73416]
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Conditional Tests: (each variable separately)
%     'F'          'p '       'R2'            'R2adj'        'Variable'
%     [ 89.486]    [0.001]    [   0.15078]    [   0.1491]    'crim'    
%     [ 75.258]    [0.001]    [   0.12992]    [  0.12819]    'zn'      
%     [ 153.95]    [0.001]    [   0.23399]    [  0.23247]    'indus'   
%     [ 15.972]    [0.001]    [  0.030716]    [ 0.028793]    'chas'    
%     [ 112.59]    [0.001]    [    0.1826]    [  0.18098]    'nox'     
%     [ 471.85]    [0.001]    [   0.48353]    [   0.4825]    'rm'      
%     [ 83.477]    [0.001]    [   0.14209]    [  0.14039]    'age'     
%     [  33.58]    [0.001]    [  0.062464]    [ 0.060604]    'dis'     
%     [ 85.914]    [0.001]    [   0.14564]    [  0.14394]    'rad'     
%     [ 141.76]    [0.001]    [   0.21953]    [  0.21798]    'tax'     
%     [ 175.11]    [0.001]    [   0.25785]    [  0.25637]    'ptratio' 
%     [ 63.054]    [0.001]    [    0.1112]    [  0.10943]    'b'       
%     [ 601.62]    [0.001]    [   0.54415]    [  0.54324]    'lstat'   
%     [  57.11]    [0.001]    [   0.10178]    [ 0.099998]    'e'       
%     [0.14346]    [ 0.71]    [0.00028455]    [-0.001699]    'n'       
% 
% 
% Used 1000 permutation of residuals under a FULL model
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'Partial F'    'p'        'Partial R2'    'Partial R2adj'    'Cum R2adj'    'Variable'
%     [   601.62]    [0.001]    [   0.54415]    [      0.54324]    [  0.54324]    'lstat'   
%     [   131.39]    [0.001]    [  0.094415]    [     0.092619]    [  0.63712]    'rm'      
%     [   62.579]    [0.001]    [  0.040063]    [     0.038158]    [   0.6767]    'ptratio' 
%     [   18.901]    [0.001]    [  0.011684]    [    0.0097226]    [  0.68784]    'dis'     
%     [   30.457]    [0.001]    [  0.017782]    [     0.015833]    [  0.70517]    'nox'     
%     [   13.492]    [0.002]    [ 0.0076849]    [     0.005716]    [  0.71236]    'chas'    
%     [   11.448]    [0.001]    [ 0.0063872]    [    0.0044157]    [  0.71826]    'b'       
%     [   8.0832]    [0.002]    [ 0.0044465]    [    0.0024712]    [  0.72221]    'zn'      
%     [   4.0555]    [0.045]    [ 0.0022172]    [    0.0002375]    [   0.7239]    'crim'    
%     [   9.9656]    [0.002]    [ 0.0053517]    [    0.0033782]    [  0.72881]    'rad'     
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







