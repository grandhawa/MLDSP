% Examples for NP-MANOVA
% 
% by David L. Jones, July-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Rattlesnakes:                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Petraitis et al. (2001) used ANCOVA to analyse observations of the age and
% size of rattlesnake (Crotalus  lepidus) collected over six years in west
% Texas.
% 
% Petraitis, P. S., S. J. Beaupre, & A. E. Dunham. 2001. ANCOVA: Nonparametric
% and randomization approaches. Pages 116-133 in S. M. Scheiner & J. Gurevitch
% (eds.), Design and Analysis of Ecological Experiments. Oxford University
% Press. Data obtained from:
% http://www.oup-usa.org/sc/0195131878/chapter7.html#anchorChapterTable
% 
% The variables are:
% site = 1 (= Boquillas) or 2 (= Grapevine Hills)
% sex  = 1 (= male) or 2 (= female)
% age  = covariate
% svl  = snout-vent length (response variable)

% Import the data:
raw = f_importCSV('rattlesnake.csv');

% Parse the variables:
snk = f_parse(raw.dat,raw.txt);
clear raw;

% Test the null hypotheis that there is no difference in body size among
% sites, taking into account the covariable 'age':
f_ancova(snk.svl,snk.site,snk.age,1000,1,0);
% 
% ==================================================
%     Nonparametric (Permutation-based) ANCOVA:
% --------------------------------------------------
% Rx*Covariate: F =   6.4704  p =  0.00600 
% Rx Effect:    F =   6.0082  p =  0.01600 
% Covariate:    F = 109.0144  p =  0.00100 
% No. of permutations = 1000 
% ==================================================
% NOTE:
% A significant 'Rx*Covariate interaction' term suggests
% the slopes of the regression lines are NOT homogeneous
% 
% -> The slopes are NOT homogeneous because there is a significant
% 'Rx*Covariate' interaction term, so the differences in body size among
% sites varies with age (= covariable).
% 
% Repeat the analysis in order to obtain a corrected form of the response
% variable with the effect of the covariable removed. Set tol=0.054 as a
% cutoff for significance levels:
ancova = f_ancova(snk.svl,snk.site,snk.age,1000,0,0.054);
% 
% Now that the response variable has been corrected for the covariable, use
% it in an ANOVA analysis:
f_npManova(f_dis(ancova.Ycor,'euc'),snk.site,1000,1);
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'        'F'        'p'    
%     'factor 1'    [ 1]    [149.85]    [149.85]    [6.264]    [0.016]
%     'residual'    [52]    [  1244]    [23.923]    [  NaN]    [  NaN]
%     'total'       [53]    [1393.8]    [   NaN]    [  NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication)
% 
% -> There is a significant difference in body size among sites even after
%    the effect of age has been removed.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Marine Plankton:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'plankton.mat' contains data regarding the size of marine plankton 
% grown in the laboratory under 8 different culture media. 
% 
% Variables:
% Rx   = 1 of 8 growth media
% day  = elapsed date of culture 
% len  = observed length of body size of each plankter

% Load the data:
load plankton.mat

% Check the number of observations within each treatment:
f_grpSize(Rx)
% 
% ans =
% 
%    174
%    173
%    175
%    165
%    174
%    175
%    175
%    171
% 
% -> this is almost a balanced design, but not quite

% Perform ANCOVA:
f_ancova(len,Rx,day,1000,1,0);
% ==================================================
%     Nonparametric (Permutation-based) ANCOVA:
% --------------------------------------------------
% Rx*Covariate: F =   2.4579  p =  0.02500 
% Rx Effect:    F =   3.1787  p =  0.00100 
% Covariate:    F = 188.1685  p =  0.00100 
% No. of permutations = 1000 
% ==================================================
% NOTE:
% A significant 'Rx*Covariate interaction' term suggests
% the slopes of the regression lines are NOT homogeneous

% Determine which growth media resulted in different growth rates:
f_ancovaPW(len,Rx,day,1000,1,'holm');
% 
% Pair-wise comparisons of SLOPES between each treatment level:
% ----------------------------------------------------------
% 1 vs. 2: F = 0.3491  p = 1.0000
% 1 vs. 3: F = 0.6199  p = 1.0000
% 1 vs. 4: F = 4.8743  p = 0.6300
% 1 vs. 5: F = 1.1149  p = 1.0000
% 1 vs. 6: F = 0.0022  p = 1.0000
% 1 vs. 7: F = 3.6920  p = 1.0000
% 1 vs. 8: F = 1.6849  p = 1.0000
% 2 vs. 3: F = 2.2006  p = 1.0000
% 2 vs. 4: F = 8.9565  p = 0.1000
% 2 vs. 5: F = 3.1329  p = 1.0000
% 2 vs. 6: F = 0.3827  p = 1.0000
% 2 vs. 7: F = 1.9523  p = 1.0000
% 2 vs. 8: F = 1.2187  p = 1.0000
% 3 vs. 4: F = 2.5021  p = 1.0000
% 3 vs. 5: F = 0.0878  p = 1.0000
% 3 vs. 6: F = 0.9018  p = 1.0000
% 3 vs. 7: F = 8.7670  p = 0.0780
% 3 vs. 8: F = 2.4804  p = 1.0000
% 4 vs. 5: F = 1.6645  p = 1.0000
% 4 vs. 6: F = 6.6794  p = 0.3910
% 4 vs. 7: F = 20.1582 p = 0.0280 *
% 4 vs. 8: F = 3.9711  p = 0.3600
% 5 vs. 6: F = 1.5826  p = 1.0000
% 5 vs. 7: F = 10.5929 p = 0.0540 *
% 5 vs. 8: F = 2.7387  p = 1.0000
% 6 vs. 7: F = 4.5239  p = 0.5060 *
% 6 vs. 8: F = 1.7048  p = 1.0000
% 7 vs. 8: F = 0.4524  p = 1.0000
% 
% Correction method: holm
% 
% -> Using the Holmes method of correcting p-values for multiple comparison
% tests, only 4 vs. 7, 5 vs. 7, and 6 vs. 7 resulted in signficanlty different
% regression slopes.
