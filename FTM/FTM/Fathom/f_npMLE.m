function result = f_npMLE(Y,X,grp,sm)
% - distance-based, maximum likelihood estimator
%
% USAGE: result = f_npMLE(Y,X,grp,sm)
% 
% Y    = matrix of UNKNOWN data (with variables corresponding to those in X)
% X    = matrix of TRAINING data (rows = obs, cols = variables)
% grp  = column vector of integers specifying group membership for rows of X
% sm   = use spatial median instead of centroid                    (default = 1)
%
% result  = structure with the following fields:
%  .grpY  = integer specifying predicted group membership of Y
%  .L     = likelihoods
%  .PP    = posterior probabilities
%  .theta = mixing proportions
%  .type  = 'spatial median' or 'centroid'
%
% SEE ALSO: f_mle, f_npDisp

% -----References:-----
% Anderson, M. J. 2006. Distance-based tests for homogeneity of multivariate
%   dispersions. Biometrics 62: 245-253
% Anderson, M. J., K. E. Ellingsen, and B. H. McArdle. 2006. Multivariate
%   dispersion as a measure of beta diversity. Ecology Letters 9(6): 683-693

% -----Author:-----
% by David L. Jones,<djones@marine.usf.edu> Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: edited to work with new version of f_pcoa

% -----Check input & set defaults:-----
if (nargin < 4), sm = 1; end % default use spatial median

rX  = size(X,1); % # obs in training set
grp = grp(:);    % force column vector

if rX ~= numel(grp), error('X & GRP need same # of rows'); end

if (size(X,2) ~= size(Y,2))
   error('X & Y must have same # columns!');
end
% -------------------------------------

% Center data:
X = X - repmat(mean(X),size(X,1),1); % center X on its column mean
Y = Y - repmat(mean(X),size(Y,1),1); % center Y based on X's mean

% Set up groups:
uGrp  = f_unique(grp);         % unique grps, unsorted
nGrp  = length(uGrp);          % number of unique groups
theta = repmat(1/nGrp,1,nGrp); % non-informative PRIORS

% -----PCoA:-----
% Combine data, get Principal Coordinates (scaled), keep negative eigenvalues:
pcoa = f_pcoa(f_dis([X;Y],'euc'),0,1,1);

% Extract eigenvectors separately for TRAINING (= U) & UNKNOWN (= V) data:
U       = pcoa.scores(1:rX,:);
V       = pcoa.scores(rX+1:end,:);
[rU,cU] = size(U);
rV      = size(V,1);

% Check for negative eigenvalues:
if (any(pcoa.evals)<0)
   error('There are negative eigenvalues!')
end

% Get Centroids of U (TRAINING SET):
C = nan(rU,cU,nGrp); % each page a copy of that group's centroid
for i = 1:nGrp
   idx = find(grp==uGrp(i)); % get indices of rows to extract
   
   % Get spatial median/centroids of each group:
   if (sm>0)
      C(:,:,i) = repmat(median(U(idx,:)),rU,1); % REAL space
   else
      C(:,:,i) = repmat(  mean(U(idx,:)),rU,1); % REAL space
   end
end

% For U: AVERAGE distance along each AXIS (rows) to each CENTROID (cols):
U  = repmat(U,[1 1 nGrp]);              % replicate into 3 pages
zU = squeeze(mean(sqrt((C - U).^2),1)); % distance, averaged across obs

% For V: ALL distances along each AXIS (rows) to each of U's CENTROIDS (cols):
V  = repmat(V,[1 1 nGrp]);                                % replicate into 3 pages
C  = permute(repmat(squeeze(C(1,:,:)),[1 1 rV]),[3 1 2]); % new C matching V's dimensions
zV = sqrt((C - V).^2);                                    % distance for each obs

% Replicate zU, so its dimensions match zV:
zU = permute(repmat(zU,[1 1 rV]),[3 1 2]);

% -----Note:-----
% ROW  = obs
% COL  = axis
% PAGE = group
% ---------------

% 'Relative' distance (along each AXIS) that accounts for each group's
% central tendency AND dispersion:
% D = abs(zV - zU);
D = abs(zV);

% Scale 'relative' distance (along each AXIS) to sum = 1 across groups,
% convert to probability:
XIrJ = 1 - (D./repmat(sum(D,3),[1 1 3])); 

% Summarize probabilities across all axes:
XIrJ = squeeze(prod(XIrJ,2));

% Likelihood each obs is from each group (using this theta):
L = XIrJ .* repmat(theta,rV,1);

% Posterior probability % scaled likelihoods to sum = 1 across groups):
PP = L ./ repmat(sum(L,2),1,nGrp);

% Update theta:
theta = mean(PP);

% Classify unknowns:
[null,grpY] = max(L,[],2);
% Make sure groups match original labels:
for i = 1:nGrp
   idx = find(grpY==uGrp(i));
   grpY(idx) = uGrp(i);
end


% Wrap results up into a structure:
result.theta = theta; % mixing proportions
result.L     = L;     % likelihoods
result.PP    = PP;    % posterior probabilities
result.grpY  = grpY;  % integers specifying predicted group membership
if (sm>0)
   result.type = 'spatial median';
else
   result.type = 'centroid';
end
% result.raw = raw(:)';
% result.con = con(:)';

