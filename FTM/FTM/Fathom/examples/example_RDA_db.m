% Example of db-RDA (Distance-based Redundancy Analysis)
% 
% by David L. Jones, Mar-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/legendre_113.mat'

% -----NOTES:-----
% This file loosely reproduces the RDA analysis of coral reef fish abundances in
% Table 11.3 of Legendre & Legendre (1998), but employs a Bray-Curtis
% distance-based approach and focuses only on the qualitative explanatory
% variable 'substrate type' and not depth.

% -----Variables:-----
% y     = abundance of 6 species of fishes from 10 sites
% y_txt = corresponding column labels
% x     = column vector of integers specifying substrate type associated with
%        each site (1=coral, 2=sand, 3=other)
% x_txt = cell array of corresponding category labels
% depth = water depth associated with each site

% -----REFERENCES:-----
% Anderson, M. J. 2001. A new method for non-parametric multivariate
%   analysis of variance. Austral Ecology 26: 32-46.
% Legendre, P., and M. J. Anderson. 1999. Distance-based redundancy analysis:
%   testing multispecies responses in multifactorial ecological experiments.
%   Ecol. Monogr. 69: 1?24.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.

% -----ANALYSIS:-----
% Load data file:
load legendre_113.mat

% Square-root transform the data to reduce the effect of more abundance species:
y_2 = f_normal(y,'2');

% Bray-Curtis dissimilarities among samples:
disBC = f_dis(y_2,'bc');

% Dummy code the substrate type (trim last column to avoid a singular matrix):
x_dum = f_dummy(x,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Use method of Legendre & Anderson, 1999:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use PCoA to obtain a Euclidean embedding of a (semi-metric) dissimilarity
% matrix, creating a response variable (U) appropriate for Euclidean
% distance-based analyses (e.g., RDA, MANOVA, etc.)
U = f_pcoa(disBC,0,1,0); % scale eigenvectors, discard neg eigenvalues
 
% Perform RDA using scaled eigenvectors + dummy codes:
rdaBC = f_rda(U.scores,x_dum,0,1000,1);
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 4.0863    p    =  0.00100 
% R2 = 0.5386   R2adj =  0.40683 
% No. of permutations = 1000 
% --------------------------------------------------

% RDA distance biplot with vectors for X,Y:
% Replace trimmed X with untrimmed version (standardize if measured on different scales):
rdaBC.X = f_stnd(f_dummy(x,0)); 
% 
H = f_rdaPlot(rdaBC,y_2,0,[0.25],0,y_txt,x_txt);
% Customize plot:
axisVar = axis;                       % get current axes bounds
axis([-1.5 1.2 axisVar(3) axisVar(4)]); % adjust x-axis

% -> Note the environmental vectors point to those sites associated with each of
%    the 3 substrate types. For the species vectors, the length of each is
%    proportional to the importance of each species in defining the ordination
%    space, and the direction depictes the gradient in each species abundance.
%    For example, all six species tend to avoid 'sand', species 5 & 6 tend to
%    prefer 'other', 3 & 4 tend to prefer 'coral', and 1 & 2 are more/less mixed
%    across 'coral/other'. Species 1 & 2 are perhaps the most important taxa as
%    their biplot vectors show the highest correlations along Axis I, which
%    accounts for 34.30% percent of the total variance (vs. 19.57% for Axis II).

% -----Sort species according to their importance:-----
% Get length of each vector:
eDis = f_dis([0 0;H.biplotY],'euc');
eDis = eDis(2:end,1);

% Show sorted list:
[null,idxVec] = sortrows(eDis,-1);
y_txt(idxVec)
% ans = 
% 
%     'sp3'
%     'sp4'
%     'sp5'
%     'sp1'
%     'sp2'
%     'sp6'
% 
% -> note how this sorted list of species corresponds with the magnitude of
%    and direction of species biplot vectors in the RDA plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Use method of Anderson, 2001:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = f_npManova(disBC,x,1000,1);
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'         'F'         'p'    
%     'factor 1'    [ 2]    [1.2553]    [0.62766]    [4.2018]    [0.001]
%     'residual'    [ 7]    [1.0457]    [0.14938]    [   NaN]    [  NaN]
%     'total'       [ 9]    [ 2.301]    [    NaN]    [   NaN]    [  NaN]
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication) 

% -> In contrast to Legendre & Anderson's (1999) method of db-RDA, Anderson's
% (2001) method of NP-MANOVA provides a correct measure of the Total SS
% resulting in a proper F-ratio. Note that use of the former method results in 
% decomposition of a semi-metric dissimilarity metric (Bray-Curtis) in order to
% calculate Sums-of-Squares. Since negative eigenvalues are involved, the Total
% SS are inflated, which affects the resulting F-ratio (because the SS
% associated with the negative eigenvalues needs to be subtracted from the SS
% associated with the positive eigenvalues). Simply ignoring those eigenvectors
% associated with negative eigenvalues, or using one of the 'correction' methods
% proposed by Legendre & Anderson (1999) still results in inflated Total SS.
% 
% However, the NP-MANOVA procedure, as described by Anderson (2001), stops short
% of providing eigenvectors to produce an ordination diagram with which to (1)
% visualize the multivariate relationships between the response and predictor
% variables and (2) assess the relative contribution of each response variable
% to the overall response-predictor relationship.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Use method of McArdle & Anderson, 2001:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_DB = f_rdaDB(disBC,size(y_2,2),x_dum,0,1000,1);
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 4.2018    p    =  0.00100 
% R2 = 0.5456   R2adj =  0.41572 
% No. of permutations = 1000 
% --------------------------------------------------

% -> This provides a method of db-RDA which directly and correctly estimates
% Total SS the same way NP-MANOVA does, side-steps issues related to negative
% eigenvalues, and provides eigenvectors for creating ordination diagrams.

% RDA distance biplot with vectors for X & Y scaled separately:
result_DB.X = f_dummy(x,0); % replace trimmed X with untrimmed version
f_rdaPlot(result_DB,y_2,0,[1 1/15],0,y_txt,x_txt);
axis([-0.5 .75 -0.75 1.05]); % customize axes

% Re-calculate biplot vectors with NO scaling:
H = f_rdaPlot(result_DB,y_2,0,[1],0,y_txt,x_txt);
% 
% Sort species according to their importance:
[null,idx] = sortrows(abs(H.biplotY),[-1 -2]);
y_txt(idx)     % show sorted list
% ans = 
%     'sp1'
%     'sp2'
%     'sp3'
%     'sp4'
%     'sp5'
%     'sp6'
H.biplotY(idx,:) % show corresponding (unscaled) values
% ans =
%       -5.1755     0.096578
%       -4.6255     -0.55542
%       -3.5383       14.982
%       -3.0584        12.95
%       -2.2808      -9.5661
%       -1.6128      -6.7643
% 
% -> note how this sorted list of species corresponds with the magnitude of
%    and direction of species biplot vectors in the RDA plot; use unscaled
%    biplot coordinates when reporting relative importance of species in, say,
%    publications, etc.
