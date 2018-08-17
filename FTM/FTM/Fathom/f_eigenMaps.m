function mem = f_eigenMaps(N,A,iter,neg)
% - create Moran's eigenvector maps from spatial coordinates
%
% USAGE: mem = f_eigenMaps(N,A,{iter},{neg}})
%
% N    = neighbor graph created by f_delaunay, f_dnn, f_gabriel, f_mst, or
%        f_relNeigh
% A    = weighting matrix (i.e., geographic similarity) created by f_dis2sim
% iter = # iterations for permutation test   (default = 0)
% neg  = keep negative eigenvalues           (default = 0)
%
% mem = structure with the following fields:
%  .evects  = positive eigenvectors
%  .evals   = positive eigenvalues
%  .MC      = Moran's coefficient of spatial autocorrelation
%  .p       = permutation-bases significance of each MC
%  .W       = spatial weighting (connectivity) matrix
%
% SEE ALSO: f_eigenMapsStepwise, f_pcnm

% -----Notes:-----
% This function creates eigenvector maps from a 2-d matrix of spatial
% coordinates based on Moran's I, following Dray et al. (2006). A spatial
% weighting matrix 'W' is created from the element-wise product of a binary
% neighborhood (connectivity) graph 'B' by a weighting matrix 'A'. Eigenvalue
% decomposition of 'W' produces a spectral decomposition of space, allowing the
% modelling of spatial structures at scales as large as the extent of the
% sampling design and as small as the largest connections in 'W' (Borcard  et al.,
% 2004). The resulting coordinates (eigenvectors) can serve as explanatory
% variables in subsequent analyses using multiple regression and canonical
% analysis (e.g., RDA). Spatial coordinates are usually in UTM (or similar)
% format rather than polar coordinates (i.e., Lat/Long).
%
% Dray et al. (2006) suggest NOT scaling the resulting eigenvectors because (1)
% those associated with negative eigenvalues will become imaginary and you will
% not be able to model negative autocorrelation to describe local structures and
% (2) they're standardized/re-scaled anyway when used in regression/canonical
% analysis, which has no effect on R^2 or the fitted values. This function
% follows their suggestion and does NOT scale the eigenvectors associated with
% eigenvector maps (as per  Dray's implementation in 'spacemakeR for R').
% 
% Note that positve eigenvalues are associated with positive autocorrelation and
% negative eigenvaleus are associated with negative autocorrelation.

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

% May-2013: updated documentation

% -----Check input & set defaults:-----
if (nargin < 2), iter = 0; end; % no permutation test by default
if (nargin < 3), neg  = 0; end; % default don't keep negative eigenvalues

if isnan(sum(N.dat(:)))
   error('N.dat is empty, re-create with original coordinates!')
end

if (size(N.B,1)~=size(A,2))
   error('Size mismatch b/n N.B and A!')
end

% Set tolerance for eigenvalues > 0:
tol = sqrt(eps);
% ----------------------

n = size(N.dat,1);     % # sites

% Spatial weighting matrix:
W = N.B .* A; % 'binary connectivity matrix' x 'weighting matrix'

% ----- PCoA modified for similarity matrices:-----
uno = ones(n,1);
I   = eye(n,n);
G   = (I-(1/n)*(uno*uno'))*W*(I-(1/n)*(uno*uno')); % Gower's centered matrix

% Eigenanalysis (eigenvectors are normalized to lengths = 1):
[evects,evals] = f_eig(G);

% Discard eigenvalues = 0:
idx    = find (abs(evals)>tol);
evects = evects(:,idx);
evals  = evals(idx);

% Only need n-1 axes for n objects:
if (size(evals,1) > n-1)
   evects = evects(:,1:(n-1));
   evals  = evals(1:(n-1));
end

% Keep negative eigenvalues:
if (neg<1)
   idx    = find(evals>0);
   evects = evects(:,idx);
   evals  = evals(idx);
end
% -------------------------------------------------

% -----Moran's Coefficient of Spatial Autocorrelation:-----
nV = size(evects,2); % get # eigenvectors
MC = zeros(nV,1);    % preallocate
p  = zeros(nV,1);

for i=1:nV
   [MC(i),p(i)] = f_moran(evects(:,i),W,iter);
end
% ---------------------------------------------------------

% -----Wrap results up into a structure:-----
mem.evects = evects;
mem.evals  = evals;
mem.MC     = MC;
mem.p      = p;
mem.W      = W;

