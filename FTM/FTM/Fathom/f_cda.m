function result = f_cda(x,y,method,iter,verb,sm)
% - Canonical Discriminant Analysis
%
% USAGE: result = f_cda(x,y,method,iter,verb,sm);
%
% x       = input data (rows = objects, cols = variables)
% y       = column vector of integers specifying group membership
% method  = center (=1, default) or stardardize (=2)
% iter    = # iterations for permutation test           (default = 0)
% verb    = optionally send result to display           (default = 1)
% sm      = use spatial median instead of centroid      (default = 0)
%
% result  = structure of outputs with the following fields:
%  .scores    = coordinates in canonical space (= canonical variates)
%  .U         = eigenvectors
%  .Cvects    = canonical eigenvectors (= canonical coefficients)
%  .evals     = canonical eigenvalues
%  .centroids = group centers in canonical space
%  .amongVar  = fraction of variance (among groups) explained
%  .lambda    = Wilks' lambda  (0 = maximum group dispersion, 1 = none)
%  .trc       = canonical trace statistic and randomized probability
%  .grs       = greatest root stat and randomized probability
%  .y         = group indicator (for f_cdaPlot)
%  .gLabels   = cell array of group labels (for f_cdaPlot)
%  .method    = type function (for f_cdaPlot)
%  .type      = 'centroid' or  'spatial 'median'
%
% SEE ALSO: f_cdaPlot, f_cdaCV

% -----Notes:-----
% This function is used to perform a classical Canonical Discriminant Analysis.
%
% METHOD = 1: variables comprising matrix X are centered, so columns of
% the canonical eigenvectors (Cvects) are CLASSIFICATION FUNCTIONS. Use this
% method when you plan to place new objects in the canonical space.
%
% METHOD = 2: variables comprising matrix X are standardized, so columns of
% the canonical eigenvectors (Cvects) are DISCRIMINANT FUNCTIONS. Use this
% method to assess the relative importance of the original variables of X
% in discriminating groups.
%
% SCORES are the coordinates of the original (centered, or standardized)
% data projected in new canonical space. They are obtained by multiplying
% the original data by the Canonical Eigenvectors.
%
% CENTROIDS are the coordinates of the group means (or medians) projected
% in the new canonical space.
%
% CVECTS (Canonical Eigenvectors or canonical coefficients) are the
% normalized orthogonal eigenvectors defining the canonical space of the
% discriminant analysis.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, Sep-2003
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% Sep-2003: uses f_eig instead of svd
% Jan-2004: tune up; added Wilks' lambda, trc, grs, permutation test, plotting
%           now handled by f_cdaPlot
% Mar-2008: updated documentation; commented out unused variables
% Feb-2011: preserved order of observations in data by removing initial
%           sorting of data; replaced 'unique' with 'f_unique'; now
%           supports spatial median; replaced find with logical indices

% -----Check input & set defaults:-----
if (nargin < 3), method =  1; end % center X by default
if (nargin < 4), iter   =  0; end % no permutation test by default
if (nargin < 5), verb   =  1; end % send output to display by default
if (nargin < 6), sm     =  0; end % default use centroid

if ((method ~= 1) && (method ~= 2)),
   error('Unknow method specified!'); end

if (size(x,1) ~= size(y,1)), error('Size mismatch b/n X & Y'); end;

% Centroid method:
if (sm>0)
   cMethod = 2; % spatial median
else
   cMethod = 1; % centroid
end
% -------------------------------------

[nRows,nCols] = size(x);

% % Sort data by groups:
% temp = sortrows([x y],nCols+1);
% x    = temp(:,1:nCols);
% y    = temp(:,nCols+1);
% clear temp;

gLabels = cellstr(num2str(f_unique(y))); % create group labels for f_cdaPlot
y       = f_recode(y);                   % recode groups to start at 1
grps    = f_unique(y);                   % unique groups, unsorted
noGrps  = size(grps(:),1);               % # of groups


if (method==1) % ---Center data on column mean---
   Xctr = f_center(x);
   
   % Centroids:
   centroids = f_centroid(Xctr,y);
   
   % Total dispersion:
   T = Xctr'*Xctr;
   
   % Within-group dispersion, pooled:
   W = zeros(nCols,nCols); % initialize
   for i=1:noGrps
      idx = (y==grps(i));
      W   = W + (f_center(x(idx,:))'*f_center(x(idx,:)));
   end
else % (method==2) % ---Standardize data on columns---
   Xstd = f_stnd(x);
   
   % Centroids:
   centroids = f_centroid(Xstd,y,cMethod);
   
   % Total dispersion:
   T = Xstd'*Xstd;
   
   % Within-group dispersion, pooled:
   W = zeros(nCols,nCols); % initialize
   for i=1:noGrps
      idx = (y==grps(i));
      W = W + (f_center(Xstd(idx,:))'*f_center(Xstd(idx,:)));
   end
end

% Divide by degrees of freedom:
% S = T/(nRows-1);
V = W/(nRows-noGrps);

% Among group dispersion:
B = T - W;
A = B/(noGrps-1);

% Wilks' lambda:
lambda = det(W)/det(T); % (eq. 11.35)

% Eigenanalyis:
[U,evals] = f_eig(f_inv(V)*A); % ratio of Among-to-Within

% Trim unnecessary axes:
if (noGrps==2)
   s = 2;
else
   s = min([nCols (noGrps-1) (nRows-1)]); % # non-zero canonical eigenvalues
end
evals = evals(1:s);
U     = U(:,1:s);

% Normalize the eigenvectors:
Cvects = U*(U'*V*U)^(-0.5);

% Canonical trace statistic:
trc.stat = sum(evals);

% Greatest root statistic:
grs.stat = evals(1);

% Project data in canonical space:
if (method<2)
   scores = Xctr*Cvects;
else
   scores = Xstd*Cvects;
end

centroids = centroids*Cvects;

% Proportion of AMONG-GROUP variation explained by each axis:
amongVar = evals./sum(evals);


%-----Permutation tests:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   randTrc = zeros(iter-1,1); % preallocate results array
   randGrs = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      % Permute grouping vector:
      randResult = f_cda(x,f_shuffle(y),method,0,0);
      randTrc(i) = randResult.trc.stat; % permuted trace stat
      randGrs(i) = randResult.grs.stat; % permuted greatest root stat
   end
   
   j1 = find(randTrc >= trc.stat); % get randomized stats >= to observed statistic
   j2 = find(randGrs >= grs.stat);
   
   trc.p = (length(j1)+1)./(iter); % count values & convert to probability
   grs.p = (length(j2)+1)./(iter);
else
   trc.p = NaN;
   grs.p = NaN;
end
%-----------------------------


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('CANONICAL DISCRIMINANT ANALYSIS:\n');
   if (method<2)
      fprintf(' Identification Functions (method = 1) \n');
   else
      fprintf(' Discriminant Functions (method = 2) \n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf('Trace stat          = %-3.4f  p =  %3.5f \n',trc.stat,trc.p);
   fprintf('Greatest root stat  = %-3.4f  p =  %3.5f \n',grs.stat,grs.p);
   fprintf('No. of permutations = %d \n',iter);
   
   fprintf('\nWilks'' lambda       = %-3.4f \n',lambda);
   fprintf('--------------------------------------------------\n');
   
   fprintf('\nCanonical Eigenvalues:\n');
   fprintf('  %-3.4f',evals);
   
   fprintf('\n\nFraction of AMONG-GROUP variance explained:\n');
   fprintf('------------------------------\n');
   fprintf('Canonical axes: \n');
   fprintf('  %-3.4f',amongVar);
   fprintf('\n  cumulative:\n');
   fprintf('  %-3.4f',cumsum(amongVar));
   
   fprintf('\n==================================================\n');
end;
% ---------------------------------


% -----Wrap results up into a structure:-----
result.scores     = scores;
result.U          = U;
result.Cvects     = Cvects;
result.evals      = evals;
result.centroids  = centroids;
result.amongVar   = amongVar;
result.lambda     = lambda;
result.trc        = trc;
result.grs        = grs;
result.y          = y;
result.gLabels    = gLabels;
result.method     = method;
if (sm>0)
   result.type = 'spatial median';
else
   result.type = 'centroid';
end
% ------------------------------------------


% Below is untested code for determining correlations with original
% variables...F_CDAPLOT uses CVECTS instead.
% % -----Correlations with variables:-----
% nr    = size(Xstd,2);        % # of variables
% nc    = size(centroids,2);   % # of canonical dimensions
% crds  = zeros(size(y,1),nc); % initialize
%
%
% for i=1:noGrps
%    idx = find(y==grps(i));
%    crds(idx,:) = repmat(centroids(i,:),length(idx),1);
% end
%
% for i=1:nr % variables
%    for j=1:nc % dimensions
%       dave(i,j) = f_corr(Xstd(:,i),crds(:,j));
%    end
% end


