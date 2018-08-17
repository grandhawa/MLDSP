% Examples for NP-MANOVA
% 
% by David L. Jones, 2002-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% updated June-2012

% 1) Example of a One-Way Model I Anova with replication (balanced design). This
% is data from Box 9.5 of Sokal & Rohlf (1995) and has 1 response variable
% (age), and 1 fixed factor (clone).
% 
% Load the data file:
load sr_09p5.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       =  f_dis(age,'euc');
result_01 = f_npManova(dis,clone,1000,1);
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'           'MS'           'F'           'p'    
%     'factor 1'    [ 1]    [0.0064286]    [0.0064286]    [0.014063]    [0.832]
%     'residual'    [12]    [   5.4857]    [  0.45714]    [     NaN]    [  NaN]
%     'total'       [13]    [   5.4921]    [      NaN]    [     NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 


% 2) Example of a Two-way Model I ANOVA with replication. This data is from
% Table 11.1 of Sokal & Rohlf (1995) and has 1 response variable (food) and 2
% fixed factors: fat and sex.
% 
% Load the data file:
load sr_11p1.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(food,'euc');
result_02 = f_npManova(dis,[sex fat],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% Select model...[O will cancel]
% 21
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'        'MS'        'F'         'p'    
%     'factor 1'      [ 1]    [3780.8]    [3780.8]    [2.5925]    [ 0.18]
%     'factor 2'      [ 1]    [ 61204]    [ 61204]    [41.969]    [0.001]
%     'factor 1x2'    [ 1]    [918.75]    [918.75]    [  0.63]    [0.471]
%     'residual'      [ 8]    [ 11667]    [1458.3]    [   NaN]    [  NaN]
%     'total'         [11]    [ 77570]    [   NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication) 


% 3) Example a Two-Way Mixed-Model ANOVA without replication. This data is from
% Box 11.3 of Sokal & Rohlf (1995) and has 1 response variable (temp) and 2
% factors: depth (fixed) and day (random).
% 
% Load the data file:
load sr_11p3.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(temp,'euc');
result_03 = f_npManova(dis,[depth day],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% Select model...[O will cancel]
% 22
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'          'F'        'p'    
%     'factor 1'    [ 9]    [2119.7]    [  235.52]    [ 2835]    [0.001]
%     'factor 2'    [ 3]    [ 0.562]    [ 0.18733]    [2.255]    [0.095]
%     'residual'    [27]    [ 2.243]    [0.083074]    [  NaN]    [  NaN]
%     'total'       [39]    [2122.5]    [     NaN]    [  NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data has NO replication) 


% 4) Example of a Two-Way Model II Nested ANOVA (balanced). This data is from
% Box 10.1 of Sokal & Rohlf (1995) and has 1 response variable (wing), and 2
% factors: cage (random) and female (nested).
% 
% Load the data file:
load sr_10p1.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(wing,'euc');
result_04 = f_npManova(dis,[cage female],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% Select model...[O will cancel]
% 23
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'        'F'         'p'    
%     'factor 1'    [ 2]    [665.68]    [332.84]    [1.7409]    [0.262]
%     'factor 2'    [ 9]    [1720.7]    [191.19]    [146.88]    [0.001]
%     'residual'    [12]    [ 15.62]    [1.3017]    [   NaN]    [  NaN]
%     'total'       [23]    [  2402]    [   NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication) 


% 5) Example of a Two-Way Model II Nested ANOVA (unbalanced). This data is from
% Box 10.6 of Sokal & Rohlf (1995) and has 1 response variable (ph) and 2
% factors: dam (random) and sire (nested).
% 
% Load the data file:
load sr_10p6.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(ph,'euc');
result_05 = f_npManova(dis,[dam sire],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% Select model...[O will cancel]
% 23
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'     'SS'        'MS'        'F'         'p'    
%     'factor 1'    [ 14]    [1780.2]    [127.16]    [3.4957]    [0.005]
%     'factor 2'    [ 22]    [800.24]    [36.374]    [1.4705]    [0.081]
%     'residual'    [123]    [3042.5]    [24.736]    [   NaN]    [  NaN]
%     'total'       [159]    [5622.9]    [   NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication)


% 6) Example of a Two-way Model III ANOVA with no replication. This data is from
% Example 12.4 of Zar (1999) and has 1 response variable (gain) and 2 random
% factors: diet and block. This is also known as a randomized block design.
% 
% Load the data file:
load zar_12p4.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(gain,'euc');
result_06 = f_npManova(dis,[diet block],1000,1);
% 
% ========================================
%   Please specify the ANOVA model 
%   for 2-way ANOVA factors 1 & 2:  
% -----------------------------------------
% 1 & 2 are fixed                    [21] 
% 1 is fixed or random, 2 is random  [22] 
% 1 is fixed or random, 2 is nested  [23] 
% Select model...[O will cancel]
% 22
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'         'F'         'p'    
%     'factor 1'    [ 3]    [27.426]    [ 9.1418]    [11.825]    [0.001]
%     'factor 2'    [ 4]    [62.647]    [ 15.662]    [20.259]    [0.001]
%     'residual'    [12]    [ 9.277]    [0.77308]    [   NaN]    [  NaN]
%     'total'       [19]    [ 99.35]    [    NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data has NO replication) 

% 7) Example of a Three-way Model I ANOVA with replication. This data is from
% Zar (1999) and has 1 response variable (rate) and 3 fixed factors: species,
% temp, and sex.
% 
% Load the data file:
load zar_14p1.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(rate,'euc');
result_07 = f_npManova(dis,[species temp sex],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% Select model...[O will cancel]
% 31
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'          'df'    'SS'           'MS'           'F'         'p'    
%     'factor 1'        [ 2]    [   1.8175]    [  0.90875]    [24.475]    [0.001]
%     'factor 2'        [ 2]    [   24.656]    [   12.328]    [332.02]    [0.001]
%     'factor 3'        [ 1]    [0.0088889]    [0.0088889]    [0.2394]    [ 0.61]
%     'factor 1x2'      [ 4]    [   1.1017]    [  0.27542]    [7.4177]    [0.001]
%     'factor 1x3'      [ 2]    [  0.37028]    [  0.18514]    [4.9863]    [0.009]
%     'factor 2x3'      [ 2]    [  0.17528]    [ 0.087639]    [2.3603]    [0.112]
%     'factor 1x2x3'    [ 4]    [  0.22056]    [ 0.055139]    [ 1.485]    [0.214]
%     'residual'        [54]    [    2.005]    [  0.03713]    [   NaN]    [  NaN]
%     'total'           [71]    [   30.355]    [      NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication) 

% 8) Example of a Three-way Model I ANOVA with no replication. This data is from
% Box 12.1 of Sokal & Rohlf (1995) and has 1 response variable (time) and 3
% fixed factors: Temp, CN, and O2
% 
% Load the data file:
load sr_12p1.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(time,'euc');
result_08 = f_npManova(dis,[temp cn o2],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% 
% Select model...[O will cancel]
% 31
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'           'MS'        'F'          'p'    
%     'factor 1'      [ 2]    [    57116]    [ 28558]    [ 441.51]    [0.001]
%     'factor 2'      [ 4]    [    55545]    [ 13886]    [ 214.68]    [0.001]
%     'factor 3'      [ 2]    [   3758.8]    [1879.4]    [ 29.055]    [0.001]
%     'factor 1x2'    [ 8]    [   3685.9]    [460.73]    [ 7.1229]    [0.002]
%     'factor 1x3'    [ 4]    [   97.067]    [24.267]    [0.37516]    [0.831]
%     'factor 2x3'    [ 8]    [   5264.5]    [658.07]    [ 10.174]    [0.001]
%     'residual'      [16]    [   1034.9]    [64.683]    [    NaN]    [  NaN]
%     'total'         [44]    [1.265e+05]    [   NaN]    [    NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data has NO replication)

% 9) Example of a Three-way Model II ANOVA with  replication. This data ships
% with MINITAB and has 1 response variable (thickness), 2 fixed factors (time,
% setting) and 1 random factor: operator.
% 
% Load the data file:
load thick.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(thickness,'euc');
result_09 = f_npManova(dis,[time setting operator],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% 
% Select model...[O will cancel]
% 32
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'          'df'    'SS'        'MS'        'F'          'p'    
%     'factor 1'        [ 1]    [     9]    [     9]    [0.29032]    [0.678]
%     'factor 2'        [ 2]    [ 15676]    [7838.2]    [ 73.178]    [0.004]
%     'factor 3'        [ 2]    [1120.9]    [560.44]    [ 4.9114]    [0.066]
%     'factor 1x2'      [ 2]    [ 114.5]    [ 57.25]    [ 2.3854]    [0.199]
%     'factor 1x3'      [ 2]    [    62]    [    31]    [ 1.2917]    [ 0.36]
%     'factor 2x3'      [ 4]    [428.44]    [107.11]    [  4.463]    [0.092]
%     'factor 1x2x3'    [ 4]    [    96]    [    24]    [  7.082]    [0.001]
%     'residual'        [18]    [    61]    [3.3889]    [    NaN]    [  NaN]
%     'total'           [35]    [ 17568]    [   NaN]    [    NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication)

% 10) Example of a Three-way Model II ANOVA with no replication. This data is
% from Brownlee (1960:p.516) and has 1 response variable (counts), 2 fixed
% factors: bottle and tube, and 1 random factor: sample.
% 
% Load the data file:
load milk.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(counts,'euc');
result_10 = f_npManova(dis,[bottle tube sample],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% Select model...[O will cancel]
% 32
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'         'MS'         'F'          'p'    
%     'factor 1'      [ 1]    [0.34722]    [0.34722]    [0.14066]    [0.705]
%     'factor 2'      [ 2]    [ 15.361]    [ 7.6806]    [ 7.6903]    [0.003]
%     'factor 3'      [11]    [ 93.486]    [ 8.4987]    [ 3.6208]    [0.075]
%     'factor 1x2'    [ 2]    [ 1.3611]    [0.68056]    [0.60767]    [0.542]
%     'factor 1x3'    [11]    [ 27.153]    [ 2.4684]    [ 2.2041]    [0.059]
%     'factor 2x3'    [22]    [ 21.972]    [0.99874]    [0.89177]    [0.604]
%     'residual'      [22]    [ 24.639]    [ 1.1199]    [    NaN]    [  NaN]
%     'total'         [71]    [ 184.32]    [    NaN]    [    NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data has NO replication) 

% 11) Example of a Three-way Model III ANOVA with replication. This data is from
% Table 23.4 of Neter et al. (1996) and has 1 response variable (tol) and 3
% random factors: gender, fat, and smoke.
% 
% Load the data file:
load exercise.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(tol,'euc');
result_11 = f_npManova(dis,[gender fat smoke],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% Select model...[O will cancel]
% 33
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'          'df'    'SS'        'MS'        'F'          'p'    
%     'factor 1'        [ 1]    [176.58]    [176.58]    [ 7.7278]    [0.052]
%     'factor 2'        [ 1]    [242.57]    [242.57]    [ 2.8797]    [ 0.12]
%     'factor 3'        [ 1]    [70.384]    [70.384]    [0.86198]    [0.252]
%     'factor 1x2'      [ 1]    [ 13.65]    [ 13.65]    [ 7.2981]    [0.217]
%     'factor 1x3'      [ 1]    [ 11.07]    [ 11.07]    [ 5.9187]    [0.268]
%     'factor 2x3'      [ 1]    [72.454]    [72.454]    [ 38.737]    [  0.1]
%     'factor 1x2x3'    [ 1]    [1.8704]    [1.8704]    [0.20036]    [0.639]
%     'residual'        [16]    [149.37]    [9.3354]    [    NaN]    [  NaN]
%     'total'           [23]    [737.95]    [   NaN]    [    NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication)


% 12) Example of a Three-way Model III ANOVA with no replication. This data is
% from Brownlee (1960:p.516) and has 1 response variable (counts) and 3 random 
% factors: bottle, tube, and sample. 
% 
% Load the data file:
load milk.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(counts,'euc');
result_12 = f_npManova(dis,[bottle tube sample],1000,1);
%  
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% Select model...[O will cancel]
% 33
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'         'MS'         'F'          'p'    
%     'factor 1'      [ 1]    [0.34722]    [0.34722]    [0.17113]    [0.546]
%     'factor 2'      [ 2]    [ 15.361]    [ 7.6806]    [ 13.731]    [0.026]
%     'factor 3'      [11]    [ 93.486]    [ 8.4987]    [ 3.6208]    [0.071]
%     'factor 1x2'    [ 2]    [ 1.3611]    [0.68056]    [0.60767]    [0.542]
%     'factor 1x3'    [11]    [ 27.153]    [ 2.4684]    [ 2.2041]    [0.055]
%     'factor 2x3'    [22]    [ 21.972]    [0.99874]    [0.89177]    [0.578]
%     'residual'      [22]    [ 24.639]    [ 1.1199]    [    NaN]    [  NaN]
%     'total'         [71]    [ 184.32]    [    NaN]    [    NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data has NO replication)

% 13) Example of a Three-way, Fully Nested, Mixed Model ANOVA (balanced). This
% data is from Box 10.5 of Sokal & Rohlf (1995) and has 1 response variable
% (gly), 1 fixed factor (rx) and 2 nested factors: rat( nested in rx) and prep
% (nested in rat).
% 
% Load the data file:
load sr_10p5.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(gly,'euc');
result_13 = f_npManova(dis,[rx rat prep],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% Select model...[O will cancel]
% 36
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'        'F'         'p'    
%     'factor 1'    [ 2]    [1557.6]    [778.78]    [ 2.929]    [0.178]
%     'factor 2'    [ 3]    [797.67]    [265.89]    [5.3715]    [0.012]
%     'factor 3'    [12]    [   594]    [  49.5]    [2.3386]    [0.048]
%     'residual'    [18]    [   381]    [21.167]    [   NaN]    [  NaN]
%     'total'       [35]    [3330.2]    [   NaN]    [   NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication)


% 14) Example of a Three-way, Cross-Nested ANOVA with replication. This data is
% from Brownlee (1960:p.530) and has 1 response variable (qual), 2 fixed
% factors: ann and loc, and 1 nested factor: coil (nested in ann). This design
% is also known as partially nested or partially hierarchical.
% 
% Load the data file:
load steel.mat
% 
% Create a Euclidean distance matrix, run the ANOVA:
dis       = f_dis(qual,'euc');
result_14 = f_npManova(dis,[ann loc coil],1000,1);
% 
% ===========================================================
%   Please specify the ANOVA model    
%   for 3-way ANOVA factors 1, 2, & 3: 
% -----------------------------------------------------------
% All factors fixed                                     [31] 
% 1 & 2 are fixed, 3 is random                          [32] 
% 1 is fixed or random, 2 & 3 are random                [33] 
% 1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] 
% 1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] 
% 3 nested in 2 nested in 1             (Fully Nested)  [36] 
% Select model...[O will cancel]
% 34
% Permuting the data 999 times...
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'        'df'    'SS'        'MS'        'F'          'p'    
%     'factor 1'      [ 1]    [  2646]    [  2646]    [  1.091]    [0.361]
%     'factor 2'      [ 1]    [1872.7]    [1872.7]    [ 35.389]    [0.001]
%     'factor 3'      [ 4]    [9701.3]    [2425.3]    [ 45.833]    [0.001]
%     'factor 1x2'    [ 1]    [16.667]    [16.667]    [0.31496]    [0.599]
%     'factor 2x3'    [ 4]    [211.67]    [52.917]    [0.50039]    [0.775]
%     'residual'      [12]    [  1269]    [105.75]    [    NaN]    [  NaN]
%     'total'         [23]    [ 15717]    [   NaN]    [    NaN]    [  NaN]
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% (Data with replication)

