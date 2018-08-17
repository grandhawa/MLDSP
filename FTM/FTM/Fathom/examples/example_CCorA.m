% Example of Canonical Correlation Analysis (CCorA)
% 
% by David L. Jones, Mar-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/doubs.mat'
% These data include the abundances of 27 species of fishes from the Doubs
% river near the border of Switzerland and France collected by Verneaux
% (1973) and included as an example data set from Borcard et al. (2011).
% This example follows that presented in their Chapter 6.
% 
% Borcard, D., F. Gillet, and P. Legendre. 2011. Numerical Ecology with R.
% Springer, NY.

% -----Variables:-----
% bio      = biotic data (30 obs, 27 species)
% bio_txt  = corresponding species names
% env      = environmental data (30 obs, 11 variables)
% env_txt  = corresponding variable labels
% site     = site number
% site_txt = corresponding label
% spa      = spatial coordinates
% spa_txt  = corresponding variable labels

% Load data:
load doubs.mat

% Remove empty site 8:
bio(8,:)    = [];
env(8,:)    = [];
spa(8,:)    = [];
site(8)     = [];
site_txt(8) = [];

% Remove 'das' from the env data:
env(:,1)   = [];
env_txt(1) = [];

% Create two subsets of env data:
envtopo     = env(:,1:3);   % physiography (upstream-downstream gradient)
envtopo_txt = env_txt(1:3);
envchem     = env(:,4:10);  % water quality
envchem_txt = env_txt(4:10);

% Transform data:
envtopo2          = envtopo; % make a copy
envtopo2(:,1:2)   = log(envtopo2(:,1:2));
envtopo2(:,3)     = f_normal(envtopo2(:,3),'2');
envchem2          = envchem; % make a copy
envchem2(:,[3 7]) = log(envchem2(:,[3 7]));
envchem2(:,4)     = f_normal(envchem2(:,4),'2');
envchem2(:,5)     = f_normal(envchem2(:,5),'ln');

% Perform CCorA using standardized variables:
result = f_CCorA(f_stnd(envchem2),f_stnd(envtopo2),1000,1);
% 
% Permuting the data 999 times...
% 
% ==================================================
%       Canonical Correlation Analysis (CCorA:
% --------------------------------------------------
% Trace Stat    = 1.8211  p =  0.00100 
% Greatest Root = 0.8482  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% Canonical Correlations:
%   0.9210  0.7637  0.6242
% Squared Canonical Correlations (= delta^2):
%   0.8482  0.5833  0.3896
% --------------------------------------------------
% Redundancy Statistics:
% Y|X: R2 = 0.4917   R2adj =  0.43075 
% X|Y: R2 = 0.7821   R2adj =  0.70946 
% ==================================================

% Create plot:
f_CCorAplot(result,f_stnd(envchem2),envchem_txt,f_stnd(envtopo2),envtopo_txt,...
   site_txt,0.1)



% -----Compare with output from R:-----
% > # CCorA (on standardized variables)
% > 
% > chem.topo.ccora <- CCorA(envchem2, envtopo2, stand.Y=TRUE, stand.X=TRUE, 
% +   nperm=999)
% > chem.topo.ccora
% 
% Canonical Correlation Analysis
% 
% Call:
% CCorA(Y = envchem2, X = envtopo2, stand.Y = TRUE, stand.X = TRUE,      nperm = 999) 
% 
%              Y X
% Matrix Ranks 7 3
% 
% Pillai's trace:  1.821107 
% 
% Significance of Pillai's trace:
% based on 999 permutations: 0.001 
% from F-distribution:  1.1352e-06 
% 
%                        CanAxis1 CanAxis2 CanAxis3
% Canonical Correlations  0.92100  0.76372   0.6242
% 
%                      Y | X  X | Y
% RDA R squares      0.49174 0.7821
% adj. RDA R squares 0.43075 0.7095
