% Example of Multivariate Correlogram
% 
% by David L. Jones, Mar-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/mites.mat'

% These data include the abundances of 35 morpho-species of orbatid mites
% from 70 cores taken from the University of Montreal's biological station
% in Quebec, Canada in June 1989. These data are included as example data
% from Borcard et al. (2011), and these examples follow that presented in
% their Chapter 7.
% 
% Borcard, D., F. Gillet, and P. Legendre. 2011. Numerical Ecology with R.
% Springer, NY.

% -----Variables:-----
% bio      = biotic data (70 obs, 35 species)
% bio_txt  = corresponding species names
% env      = environmental data (70 obs, 5 variables):
%            substrate density (g/dm^-3)
%            water content (g/dm^-3)
%            substrate       = 1: sphagn1, 2: sphagn2, 3: sphagn3, 4: sphagn4,
%                              5: interface, 6: barepeat, 7: litter
%            shrubs          = 0: none, 1: few, 2: many
%            microtopography = 1: blanket, 2: hummock 
% env_txt  = corresponding variable labels
% site     = site number
% site_txt = corresponding label
% spa      = spatial coordinates
% spa_txt  = corresponding variable labels

% Load the data:
load mites;

% Hellinger transform the biotic data:
bio_H = f_transform(bio,'hel');

% Detrend biotic data so they are 2nd order stationary:
rda_spa = f_rda(bio_H,spa,0,1000,1);
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
bio_H_dt = rda_spa.res; % residuals serve as detrended data

% Create distance matrix:
yDis = f_dis(bio_H_dt,'euc');

% Multivariate Mantel Correlogram using Sturge's Rule and Holmes correction
% for multiple testing: 
res = f_correlogram(yDis,spa,0,1000,'holm',1,1,1,1);
% 
% ============================
%     Mantel Correlogram:
% ============================
% Class:  r:     p:     pA:   
%      1 +0.1357 0.0010 0.0010 
%      2 +0.1182 0.0010 0.0020 
%      3 +0.0378 0.0570 0.0570 
%      4 -0.0986 0.0010 0.0040 
%      5 -0.1127 0.0010 0.0050 
%      6 -0.1076 0.0010 0.0060 
%      7 -0.0223 0.1250 0.1250 
%      8 +0.0234 0.1240 0.2480 
%      9 +0.0407 0.0370 0.1480 
%     10 +0.0523 0.0230 0.1150 
%     11 +0.0277 0.1570 0.3720 
%     12 -0.0521 0.0460 0.2300 
%     13 -0.0834 0.0030 0.0240 
% 
% # iterations       = 1000
% Correction Method = holm
% Progressive       = 1
% ----------------------------
% Class = distance class (half = 4.809366e+00) 
% r     = Mantel statistic 
% p     = permutation-based p-value 
% pA    = adjusted p-value 

