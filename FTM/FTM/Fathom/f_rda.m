function result = f_rda(Y,X,W,iter,verb,perm,verbPerm)
% - Redundancy Analysis (RDA)
%
% USAGE: result = f_rda(Y,X,W,iter,verb)
%
% Y    = matrix of response variables       (rows = sites, cols = species)
%
% X    = (1) matrix of explanatory variables, or
%        (2) ANOVA design matrix specified by dummy coding
%
% W    = matrix of covariables               (0 = none)
% iter = # iterations for permutation test   (default = 0)
% verb = optionally send results to display  (default = 1)
%
% result = structure returning the following values:
%  .F          = F statistic
%  .p          = p-value for F stat
%  .R2         = Coefficient of determination (R^2)
%  .R2adj      = adjusted R^2
%  .dfR        = df regression
%  .dfE        = df error
%  .dfT        = df total
%  .SSr        = SS regression
%  .SSr_W      = SS regression of covariable
%  .SSe        = SS error
%  .SSt        = SS total
%  .MSr        = MS regression
%  .MSe        = MS error
%  .X          = standardized explanatory variables
%  .Y          = centered response variables
%  .fit        = fitted values of Y
%  .res        = residual values of Y
%  .siteScores = site scores        (= WA scores)
%  .fitScores  = fitted site scores (= LC scores)
%  .resScores  = residual scores
%  .U          = canonical eigenvectors
%  .Ures       = residual eigenvectors
%  .evals      = canonical eigenvalues
%  .evalsRes   = residual eigenvalues
%  .canVar     = fraction of variance explained by canonical axes (= R^2)
%  .canCorr    = Species-Environ correlations (r) of each canonical axis
%
% SEE ALSO: f_rdaPlot, f_rdaDB, f_rdaAnova, f_rdaAIC, f_rdaStepwise, f_npManova


% -----Notes:-----
% If Y is a column vector (= single response variable), RDA is simply multiple
% linear regression analysis and the PCA steps are skipped.

% -----Scaling:-----
% In this implementation of RDA the eigenvectors are scaled to lengths of 1,
% creating DISTANCE biplots which preserve the distances among sites. These are
% interpreted as follows:
%
%     A) Distances among sites approximate their Euclidean distance.
%     B) Projecting sites onto a Y arrow approximates the value of that
%        site along the variable.
%     C) Angles among Y variables are meaningless.
%     D) Angles b/n Y and X arrows reflect their correlations.
%     E) Distances among centroids and between centroids and sites
%        approximate Euclidean distance.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% Legendre, P. 2007. Studying beta diversity: ecological variation partitioning
%   by multiple regression and canonical analysis. Journal of Plant Ecology,
%   doi:10.1093/jpe/rtm001
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.
% Peres-Neto, P. R., P. Legendre, S. Dray, & D. Borcard. 2006. Variation
%   partitioning of species data matrices: estimation and comparison of
%   fractions. Ecology 87(10): 2614-2625.

% -----Author:-----
% by David L. Jones, Mar-2003
% 
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.
%
% This code is provided as is, with no guarantees

% Oct-2003: updated documentation, fixed 'Check scaling', added verbPerm to
%           reduce clutter when run by f_rdaStepwise
% Jan-2008: updated documentation
% Feb-2008: updated documentation; changed | to ||, changed & to &&;
%           preallocated corrX_fit, corrY_fit,canonCorr; num2Str canged to
%           num2str; moved .corrX, .corrY, and .biplotX and all plotting
%           routines to 'f_rdaPlot' (as a result the 5th command-line option has
%           been removed so adjust legacy code calling this)
% Mar-2008: overhaul: F-stat relies on SS, not permuted eigenvectors;
%           minimized calculations during permutation runs; trim residual axes
%           with 0 variance; added R2adj; added ANOVA stats; PCA is skipped for
%           univarite response data; make sure W ~= X
% May-2008: added SSe, SSt to output when perm>0
% Mar-2013: modified to be compatible with changes made to f_pca
% May-2014: now outputs scores and canVar for univariate data for plotting

% -----Set defaults & check input:-----
if (nargin < 3), W        = 0; end % no covariables by default
if (nargin < 4), iter     = 0; end % no permutation test by default
if (nargin < 5), verb     = 1; end % send output to display by default
if (nargin < 6), perm     = 0; end % internal flag: not a permutation run
if (nargin < 7), verbPerm = 1; end % internal flag: verbose permutation indicator

[n,ncY] = size(Y);   % # sites, # species
ncX     = size(X,2); % # predictors
if n   ~= size(X,1), error('Y & X need same # of rows'); end;

% Check for covariables:
[nrW,ncW] = size(W);
if (nrW==1) && (ncW==1) && (sum(W(:))==0)
   partial = 0; % no partial analysis
   ncW     = 0; % # covariables
else
   partial = 1; % do partial analysis
   if (nrW ~= n)
      error('W & X need same # of rows')
   end
   if (ncW==ncX)
      if (sum(sum(X==W)) == n)
         error('X & W cannot be the same matrix!')
      end
   end
end
% ----------------------

Y = f_center(Y); % center response variables           (L&L,1998: p.580)
X = f_stnd(X);   % standardize explanatory variables
W = f_stnd(W);   % standardize covariables

% =========================================================================
if (partial<1) % No covariables:
   B    = X\Y;      % regression coefficients via QR
   Yfit = X*B;      % fitted values of Y
   Yres = Y - Yfit; % residuals

else
   % Variables + covariables:
   XW      = [X W];
   B_XW    = XW\Y;
   Yfit_XW = XW*B_XW;
   Yres    = Y - Yfit_XW;

   % Covariables:
   B_W    = W\Y;
   Yfit_W = W*B_W;
   Yres_W = Y - Yfit_W;       % needed for permutation test

   Yfit   = Yfit_XW - Yfit_W; % partial out covariables
end
% =========================================================================

% McArdle & Anderson (2001):
SSr = trace(Yfit'*Yfit); % Sum-of-Squares of regression model
SSe = trace(Yres'*Yres); % Sum-of-Squares error

if (partial>0)
   SSr_W = trace(Yfit_W'*Yfit_W); % Sum-of-Squares of covariable
else
   SSr_W = NaN;
end

% L&L,1998 (eq. 11.19); L,2007 (eq. 3):
dfR = ncX;         % Degrees-of-Freedom regression model
dfE = n-ncX-ncW-1; % Degrees-of-Freedom error
MSr = SSr/dfR;     % Mean Squares regression model
MSe = SSe/dfE;     % Mean Squares error
F   = MSr/MSe;     % Observed F-stat

% Coefficient of determination:
SSt   = trace(Y'*Y);                    % Sum-of-Squares total
R2    = SSr/SSt;                        % Legendre,2007 (eq. 5)
R2adj = 1 - ((1-R2)*[(n-1)/(n-ncX-1)]); % Legendre,2007 (eq. 6)

% If this is a permutation run, stop here:
if (perm>0)
   result.F     = F;
   result.p     = NaN;
   result.R2    = R2;
   result.R2adj = R2adj;
   result.SSe   = SSe;
   result.SSt   = SSt;
   return;
end

if ncY==1 % Univariate response data:
   siteScores = Y;
   fitScores  = Yfit;
   resScores  = Yres;
   U          = NaN;
   Ures       = NaN;
   evals      = NaN;
   evalsRes   = NaN;
   canVar     = 1;
   canonCorr  = NaN;

else % Multivariate response data:
   % PCA:
   pca      = f_pca(Yfit); % fitted Y
   pcaRes   = f_pca(Yres); % residuals
   U        = pca.evects;
   evals    = pca.evals(:)';
   Ures     = pcaRes.evects;
   evalsRes = pcaRes.evals(:)';
   clear pca pcaRes;
      
   % Project in canonical space:
   siteScores = Y*U;       % site scores        (in space of Y)         (eq.11.12)
   fitScores  = Yfit*U;    % fitted site scores (in space of X)         (eq.11.13)
   resScores  = Yres*Ures; % residual scores    (in space of residuals) (fig.11.2)

   % Trim unnecessary axes:
   s          = min([ncX ncY (n-1)]); % # non-zero canonical eigenvalues
   sRes       = min([ncY (n-1)]);     % # non-zero residual eigenvalues
   siteScores = siteScores(:,1:s);
   fitScores  = fitScores(:,1:s);
   resScores  = resScores(:,1:sRes);
   U          = U(:,1:s);
   Ures       = Ures(:,1:sRes);
   evals      = evals(1:s);
   evalsRes   = evalsRes(1:sRes);

   % Only keep residual axes with variances larger than 0 (L&L,1998:p.590)
   idx           = find(var(resScores) < 0.000001 == 1);
   evalsRes(idx) = [];
   Ures(:,idx)   = [];

   % Species-Environment correlation (r) of each exis:
   canonCorr = zeros(1,s); % preallocate
   for i = 1:s
      canonCorr(i) = f_corr(siteScores(:,i),fitScores(:,i));
   end
     
   % Proportion of variance:
   totVar = sum([evals evalsRes]); % total variance in Y (= total inertia)
   canVar = evals./totVar;         % fraction of variance explained by canonical axes
   resVar = evalsRes./totVar;      % fraction of variance explained by residual axes
end


% =========================================================================
if (iter>0)
   Fperm = zeros(iter-1,1); % preallocate result array

   if (partial<1) % NO COVARIABLES:
      % Permutation of residuals under a full model (L&L,1998 p.608):
      if (verbPerm>0)
         fprintf('\nPermuting residuals under a full model %d times...\n',iter-1);
      end

      for i = 1:(iter-1) % observed value is considered a permutation
         YresPerm   = f_shuffle(Yres,4); % permute rows of residuals
         resultPerm = f_rda(YresPerm,X,0,0,1,1);
         Fperm(i)   = resultPerm.F;
      end

   else % WITH COVARIABLES:
      % Permutation of residuals under a reduced model (L&L,1998 p.609-610):
      if (verbPerm>0)
         fprintf('\nPermuting residuals under a reduced model %d times...\n',iter-1);
      end

      for i = 1:(iter-1) % observed value is considered a permutation
         Yres_Wperm = f_shuffle(Yres_W,4); % permute rows of covariable residuals
         Yperm      = Yfit_W + Yres_Wperm; % add permuted residuals back to form new Y
         resultPerm = f_rda(Yperm,X,W,0,0,1,1);
         Fperm(i)   = resultPerm.F;
      end

   end

   j = find(Fperm >= F);      % get permuted stats >= to observed statistic
   p = (length(j)+1)./(iter); % count values & convert to probability

else
   p = NaN;
end
% =========================================================================


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   if (partial<1)
      fprintf('REDUNDANCY ANALYSIS:\n');
   else
      fprintf('Partial REDUNDANCY ANALYSIS:\n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf(' F = %-3.4f    p    =  %3.5f \n',F,p);
   fprintf('R2 = %-3.4f   R2adj =  %3.5f \n',R2,R2adj);
   fprintf('No. of permutations = %d \n',iter);
   fprintf('--------------------------------------------------\n');

   if (ncY==1) % Univariate response variable:
      fprintf('Response variable (Y) is univariate. \n');
      fprintf('--------------------------------------------------\n');
   else        % Multivariate response variabel:

      fprintf('\nCanonical Eigenvalues:\n');
      fprintf('  %-3.4f',evals);
      fprintf('\nResidual Eigenvalues:\n');
      fprintf('  %-3.4f',evalsRes);

      fprintf('\n\nSpecies-Environment Correlations (r):\n');
      fprintf('  %-3.4f',canonCorr);

      if (partial<1)
         fprintf('\n\nFraction of variance explained:\n');
      else
         fprintf('\n\nFraction of RESIDUAL variance explained (co-variable removed):\n');
      end
      fprintf('------------------------------\n');
      fprintf('Canonical axes (total = %-3.4f): \n',sum(canVar));
      fprintf('  %-3.4f',canVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %-3.4f',cumsum(canVar));
      fprintf('\n\nResidual axes  (total = %-3.4f):\n',sum(resVar));
      fprintf('  %-3.4f',resVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %-3.4f',cumsum(resVar));
      fprintf('\n==================================================\n');
   end
end
% ---------------------------------

% -----Wrap results up into a structure:-----
result.F          = F;
result.p          = p;
result.R2         = R2;
result.R2adj      = R2adj;

result.dfR        = dfR;
result.dfE        = dfE;
result.dfT        = n-1;

result.SSr        = SSr;
if (partial>0)
   result.SSr_W   = SSr_W;
end
result.SSe        = SSe;
result.SSt        = SSt;

result.MSr        = MSr;
result.MSe        = MSe;

result.X          = X;
result.Y          = Y;
result.fit        = Yfit;
result.res        = Yres;

result.siteScores = siteScores;
result.fitScores  = fitScores;
result.resScores  = resScores;

result.U          = U;
result.Ures       = Ures;
result.evals      = evals;
result.evalsRes   = evalsRes;

result.canVar     = canVar;
result.canCorr    = canonCorr;
% ------------------------------------------
