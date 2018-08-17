% Example of Nonparametric Analysis of Multivariate Dispersion 
% 
% by David L. Jones, Feb-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             VARESPEC:               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FILE: varespec.mat
% The 'varespec' data set is from the VEGAN package and consists of vegetation
% data in lichen pastures from Väre et al. (1995). There are 24 observations of
% estimated cover values for 44 species. [Väre, H., Ohtonen, R. and Oksanen, J.
% (1995) Effects of reindeer grazing on understorey vegetation in dry Pinus
% sylvestris forests. Journal of Vegetation Science 6, 523-530.]

% Variables:
% Y     = species abundances (rows = obs, cols = species)
% Y_txt = corresponding species labels
% X     = grouping variable (1 = grazed, 2 = ungrazed)
% X_txt = grouping labels as a cell array

% Load the data:
load varespec.mat

% Square-root transform the data:
Y_2 = f_normal(Y,'2');

% Bray-Curtis dissimilarities among samples:
disBC = f_dis(Y_2,'bc');

% Calculate multivariate dispersions:
md = f_npDisp(disBC,X,1,1);
% 
% ==================================================
%       Quantify Multivariate Dispersions:
% ==================================================
% 
% # Pos Eigenvalues = 17
% # Neg Eigenvalues = 6
% 
% Average distance to spatial median:
% Group 1 = 0.2925
% Group 2 = 0.2319
% ---------------------------------------------------

% Test for significant differences in among-group dispersion:
f_npManova(f_dis(md.z,'euc'),X,1000,1);

% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'          'MS'           'F'         'p'    
%     'factor 1'    [ 1]    [0.019571]    [ 0.019571]    [2.2577]    [0.164]
%     'residual'    [22]    [ 0.19071]    [0.0086685]    [   NaN]    [  NaN]
%     'total'       [23]    [ 0.21028]    [      NaN]    [   NaN]    [  NaN]
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 
% 
% -> Note: there are NO significant differences in dispersion among groups.

% Create Plot showing bootstrapped 95% confidence intervals:
hdl = f_npDispPlot(md,0.95,1000);
legend(hdl,{'grazed' 'ungrazed'}); % create legend


% Test for significant differences in central tendency:
f_npManova(disBC,X,1000,1);
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'         'MS'          'F'         'p'    
%     'factor 1'    [ 1]    [0.49898]    [ 0.49898]    [5.7625]    [0.001]
%     'residual'    [22]    [  1.905]    [0.086591]    [   NaN]    [  NaN]
%     'total'       [23]    [  2.404]    [     NaN]    [   NaN]    [  NaN]
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 
% 
% -> Note: there ARE significant differences in the abundance/species
% composition of species among the two groups (i.e., grazed vs. ungrazed)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             SPARROWS:               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FILE: sparrows.mat
% 
% The sparrow data set consists of the length of body measurements (in mm) of 49
% female sparrows made by Bumpus (1898). Birds 1 to 21 survived a severe storm
% near Brown University in Rhode Island while 22-49 died. [Bumpus, H. C. 1898.
% Biological Lectures, 11th Lecture, Marine Biology Laboratory, Woods Hole, MA,
% pp.209-226.]
% 
% 
% X = body measurements of female sparrows (in mm):
%     col 1  = total length
%     col 2  = alar length 
%     col 3  = length of beak and head, 
%     col 4  = length of humerus
%     col 5  = length of keel and sternum
% 

% Load data:
load sparrows.mat;

% Standardize each variable to z-scores (normalize):
Xs = f_stnd(X);

% Calculate multivariate dispersions:
md = f_npDisp(f_dis(Xs,'euc'),Y,0,1);
% 
% ==================================================
%       Quantify Multivariate Dispersions:
% ==================================================
% 
% # Pos Eigenvalues = 5
% # Neg Eigenvalues = 0
% 
% Average distance to centroid:
% Group 1 = 1.7309
% Group 2 = 2.2309
% ---------------------------------------------------

% Test for significant differences in among-group dispersion:
f_npManova(f_dis(md.z,'euc'),Y,10000,1);
% 
% Permuting the data 9999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'         'F'         'p'     
%     'factor 1'    [ 1]    [2.9997]    [ 2.9997]    [3.8692]    [0.0515]
%     'residual'    [47]    [36.438]    [0.77528]    [   NaN]    [   NaN]
%     'total'       [48]    [39.438]    [    NaN]    [   NaN]    [   NaN]
% 
%       # iterations =       10000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 
% 
% -> Note: there ARE significant differences in dispersion among groups.

% Create Plot:
f_npDispPlot(md,0.95,1000);
hdl = f_npDispPlot(md);
legend(hdl,{'survived' 'died'})

