function result = f_pcnm(X)
% - Principal Coordinates of Neighbor Matrices (PCNM) from spatial coordinates
%
% USAGE: result = f_pcnm(X)
%
% X = 2 column matrix of site coordinates    (e.g., [e n])
%
% result = structure with the following fields:
%  .evects  = eigenvectors
%  .evals   = positive eigenvalues
%  .dis     = euclidean distance matrix
%  .tDis    = truncated distance matrix
%
% SEE ALSO: f_eigenMaps

% -----Notes:-----
% This function creates Principal Coordinates of Neighbor Matrices (PCNM) from a
% 2-d matrix of spatial coordinates using a truncated distance matrix sensu
% Borcard & Legegendre, 2002. The resulting coordinates (eigenvectors) can serve
% as explanatory variables in subsequent analyses using multiple regression and
% canonical analysis (e.g., RDA). Spatial coordinates are usually in UTM (or
% similar) format rather than polar coordinates (i.e., Lat/Long). The procedure
% involves the eigenvalue decomposition of a neighborhood (connectivity) graph,
% in this case a truncated symmetrical distance matrix of the spatial
% coordinates. This produces a spectral decomposition of space, allowing the
% modelling of spatial structures at scales as large as the extent of the
% sampling design and as small as the truncation distance (Borcard et al.,
% 2004).
%
% Dray et al. (2006) suggest NOT scaling the resulting eigenvectors because (1)
% those associated with negative eigenvalues will become imaginary and you will
% not be able to model negative autocorrelation to describe local structures and
% (2) they're standardized/re-scaled anyway when used in regression/canonical
% analysis, which has no effect on R^2 or the fitted values. This function
% follows their suggestion and does NOT scale the eigenvectors associated with
% PCNM's (as per Dray's implementation in 'spacemakeR for R'). However, only
% positve eigenvalues (associated with positive autocorrelation) are returned.
% 
% Dray et al. (2006) explore the theoretical underpinnings of PCNM's and define
% them as a special case of the more general framework of Moran's Eigenvector
% Maps (MEM's). Thus, PCNM's are MEM's constructed using a minimum spanning tree
% (MST) as the binary connectivity matrix ('B') and a weighting matrix ('A')
% based on a truncated distance matrix. Griffith and Peres-Neto (2006) describe
% PCNM's as distance-based MEM's and construct them essentially in this manner.
% Dray et al. (2006) describe a variety of other possible ways to construct
% MEM's based on different types of binary neighborhood graphs ('B'), not just
% MST's, and describe other methods of creating the weighting matrix ('A').
% These methods are available in 'f_eigenMaps' and that function should be
% preferred over 'f_pcnm', which is included mainly for making comparisons. 

% -----References:-----
% Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological
%   data by means of principal coordinates of neighbor matrices. Ecological
%   Modelling 153: 51-68.
% Borcard, D., P. Legendre, C. Avois-Jacquet, and H. Tuomisto. 2004. Dissecting
%   the spatial structure of ecological data at multiple scales. Ecology 85(7):
%   1826-1832.
% Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial-modelling: a
%   comprehensive framework for principal coordinate analysis of neighbor
%   matrices (PCNM). Ecological Modelling 196: 483-493.
% Griffith, D. A. and P. R. Peres-Neto. 2006. Spatial modeling in ecology: the
%   flexibility of eigenfunction spatial analysis. Ecology 87(10): 2603-2613.

% -----Author:-----
% by David L. Jones, Apr-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2011: replaced f_euclid with f_dis
% Jun-2012: updated call to f_pcoa

% -----Check input & set defaults:-----
if (size(X,2)~=2) % this could be modified to handle 3-D (or more) data
   error('X must be a 2 column matrix specifying x,y coordinates of sites!');
end
% ----------------------

% Minimum Spanning Tree:
mst = f_mst(f_dis(X,'euc'),X);

% PCoA (don't scale eigenvectors, don't keep neg eigenvalues):
temp = f_pcoa(mst.tDis,0,0,0);

% -----Wrap results up into a structure:-----
result.evects = temp.scores;
result.evals  = temp.evals;
result.dis    = mst.dis;
result.tDis   = mst.tDis;
