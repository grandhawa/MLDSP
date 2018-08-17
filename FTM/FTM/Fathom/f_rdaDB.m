function result = f_rdaDB(yDis,ncY,X,W,iter,verb,perm,verbPerm)
% - distance-based Redundancy Analysis (db-RDA)
%
% USAGE: result = f_rdaDB(yDis,ncY,X,W,iter,verb);
%
% yDis = square symmetric distance matrix derived from response variables
% ncY  = # columns of original transformed data used to derive yDis
%
% X    = (1) matrix of explanatory variables, or
%        (2) ANOVA design matrix specified by dummy coding
%
% W     = matrix of covariables                                       (0 = none)
% iter  = # iterations for permutation test                        (default = 0)
% verb  = optionally send results to display                       (default = 1)
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
%  .Y          = NaN
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
% SEE ALSO: f_rdaDB_AIC, f_rda, f_rdaPlot, f_npManova, f_mregress

% -----Notes:-----
% Distance-based Redundancy Analysis is a constrained form of PCoA. Therefore,
% it may be instructive to compare an ordination plot derived from PCoA, which
% depicts a maximum of the total variation in the response variable with a
% minimum number of axes, vs. db-RDA, which depicts a maximum of the total
% variation in the repsomse variable explained by the predictor variable(s) with
% a minimum number of axes.

% -----References:-----
% Anderson, M. J. 2001. A new method for non-parametric multivariate
%   analysis of variance. Austral Ecology 26: 32-46.
% Anderson, M. J. 2002. DISTML v.2: a FORTRAN computer program to calculate a
%   distance-based multivariate analysis for a linear model. Dept. of Statistics
%   University of Auckland. (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% Legendre, P. 2007. Studying beta diversity: ecological variation partitioning
%   by multiple regression and canonical analysis. Journal of Plant Ecology,
%   doi:10.1093/jpe/rtm001
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.
%
%  Most of the equations are from McArdle & Anderson (2001).

% -----Author:-----
% by David L. Jones, Mar-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2010: sub_decompose no longer trims its results, eigenvectors are now
%           scaled
% Dec-2012: 'I' and 'uno' defined for permutation runs, too.

% -----Set defaults & check input:-----
if (nargin < 4), W        = 0; end % no covariables by default
if (nargin < 5), iter     = 0; end % no permutation test by default
if (nargin < 6), verb     = 1; end % send output to display by default
if (nargin < 7), perm     = 0; end % internal flag: not a permutation run
if (nargin < 8), verbPerm = 1; end % internal flag: verbose permutation indicator

n    = size(yDis,1); % # sites
%ncY = provided on command line
ncX  = size(X,2);    % # predictors
if n ~= size(X,1), error('yDis & X need same # of rows'); end

% Check input if this isn't a permutation run:
if (perm==0)
   if (f_issymdis(yDis) == 0)
      error('Input yDIS must be a square symmetric distance matrix!');
   end
end

% Check for covariables:
[nrW,ncW] = size(W);
if (nrW==1) && (ncW==1) && (sum(W(:))==0)
   partial = 0; % no partial analysis
   ncW     = 0; % # covariables
else
   partial = 1; % do partial analysis
   if (nrW ~= n)
      error('W & X need same # of rows');
   end
   if (ncW==ncX)
      if (sum(sum(X==W)) == n)
         error('X & W cannot be the same matrix!')
      end
   end
end
% -------------------------------------

% Explanatory/Covariables:
X = f_stnd(X);   % standardize explanatory variables
W = f_stnd(W);   % standardize covariables

% Centering and Intercept terms:
I   = eye(n,n);
uno = ones(n,1);

if ~(perm ==-1)
   % Response variable:
   A = -0.5*(yDis.^2);
   G = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno')); % Gower's centered matrix
else
   % Permutation of residuals, yDis is already 'G':
   G   = yDis;
end

% =========================================================================
if (partial<1) % No covariables:
   [Q1,R1] = qr([uno X],0); H = Q1*Q1';      % Hat-matrix
   SSr     = trace(H*G*H);                   % SS regression
   SSe     = trace((I-H)*G*(I-H));           % SS error
   
else
   % Variables + covariables:
   [Q1,R1] = qr([uno X W],0); H_XW = Q1*Q1'; % Hat-matrix
   SSr_XW  = trace(H_XW*G*H_XW);             % SS regression
   SSe     = trace((I-H_XW)*G*(I-H_XW));     % SS error
   % XW      = [X W];
   % B_XW    = XW\Y;
   % Yfit_XW = XW*B_XW;
   % Yres    = Y - Yfit_XW;
   
   % Covariables:
   [Q1,R1] = qr([uno W],0); H_W = Q1*Q1';    % Hat-matrix for
   SSr_W   = trace(H_W*G*H_W);               % SS regression for W
   SSe_W   = trace((I-H_W)*G*(I-H_W));       % SS error for W
   Gfit_W  = (H_W*G*H_W);                    % fitted values for W
   Gres_W  = G - Gfit_W;                     % residuals for W
   % Gres_W = (I-H_W)*G*(I-H_W);             % alternative formulation
   % B_W    = W\Y;
   % Yfit_W = W*B_W;
   % Yres_W = Y - Yfit_W; % needed for permutation test
   
   SSr = SSr_XW - SSr_W;                     % partial out covariables
   % Yfit   = Yfit_XW - Yfit_W; % partial out covariables
   
   % Hat matrix when partial out covariables ?????
   H = H_XW - H_W; % SEEMS RIGHT, CHECK AGAINST F_RDA OUTPUT [DLJ: AUG-2009]!!!!
end
% =========================================================================

F = (SSr/ncX)/(SSe/(n-ncX-ncW-1)); % L&L,1998 (eq. 11.19); L,2007 (eq. 3)
% -> this is the same as: F = MSr/MSe

% Coefficient of determination:
SSt   = trace(G);                       % SS total
R2    = SSr/SSt;                        % Legendre,2007 (eq. 5)
R2adj = 1 - ((1-R2)*((n-1)/(n-ncX-1))); % Legendre,2007 (eq. 6)

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

% More stats:
dfR = ncX;         % Degrees-of-Freedom regression model
dfE = n-ncX-ncW-1; % Degrees-of-Freedom error
MSr = SSr/dfR;     % Mean Squares regression model
MSe = SSe/dfE;     % Mean Squares error

% Eigenvalue decomposition:
[U,evals]       = sub_decompose(H*G*H);         % fitted Y
[Ures,evalsRes] = sub_decompose((I-H)*G*(I-H)); % residuals

% Project PCoA's in canonical space:
siteScores = G*U;                  % site scores        (in space of Y) (eq.11.12)
fitScores  = (H*G*H)*U;            % fitted site scores (in space of X) (eq.11.13)
resScores  = (I-H)*G*(I-H) * Ures; % residual scores    (in space of residuals) (fig.11.2)

% Scale axes to sqrt of their eigenvalue:
siteScores = siteScores .* repmat(abs((evals.^0.5)),size(siteScores,1),1);
fitScores  = fitScores  .* repmat(abs((evals.^0.5)),size(fitScores,1),1);
resScores  = resScores  .* repmat(abs((evalsRes.^0.5)),size(resScores,1),1);

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

% =========================================================================
if (iter>0)
   Fperm = zeros(iter-1,1); % preallocate result array
   
   if (partial<1) % NO COVARIABLES:
      % Permutation of residuals under a full model (L&L,1998 p.608):
      
      if (verbPerm>0)
         fprintf('\nPermuting residuals under a full model %d times...\n',iter-1);
      end
      
      for i = 1:(iter-1) % observed value is considered a permutation
         yDisPerm   = f_shuffle(yDis,2); % permute symmetric distance matrix
         resultPerm = f_rdaDB(yDisPerm,ncY,X,0,0,0,1,1);
         Fperm(i)   = resultPerm.F;
      end
      
   else % WITH COVARIABLES:
      % Permutation of residuals under a reduced model (L&L,1998 p.609-610):
      if (verbPerm>0)
         fprintf('\nPermuting residuals under a reduced model %d times...\n',iter-1);
      end
      
      for i = 1:(iter-1) % observed value is considered a permutation
         Gperm      = Gfit_W + f_shuffle(Gres_W); % add permuted residuals back to form new G
         resultPerm = f_rdaDB(Gperm,ncY,X,W,0,0,-1,1);
         Fperm(i)   = resultPerm.F;
      end
   end
   
   j  = find(Fperm >= F);      % get permuted stats >= to observed statistic
   p  = (length(j)+1)./(iter); % count values & convert to probability
   
else
   p  = NaN;
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
   fprintf(' F = %3.4f    p    =  %3.5f \n',F,p);
   fprintf('R2 = %3.4f   R2adj =  %3.5f \n',R2,R2adj);
   fprintf('No. of permutations = %d \n',iter);
   fprintf('--------------------------------------------------\n');
   
   if (ncY==-1) % Univariate response variable:
      fprintf('Response variable (Y) is univariate. \n');
      fprintf('--------------------------------------------------\n');
   else        % Multivariate response variabel:
      
      fprintf('\nCanonical Eigenvalues:\n');
      fprintf('  %3.4f',evals);
      fprintf('\nResidual Eigenvalues:\n');
      fprintf('  %3.4f',evalsRes);
      
      fprintf('\n\nSpecies-Environment Correlations (r):\n');
      fprintf('  %3.4f',canonCorr);
      
      if (partial<1)
         fprintf('\n\nFraction of variance explained:\n');
      else
         fprintf('\n\nFraction of RESIDUAL variance explained (co-variable removed):\n');
      end
      fprintf('------------------------------\n');
      fprintf('Canonical axes (total = %-3.4f): \n',sum(canVar));
      fprintf('  %3.4f',canVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %3.4f',cumsum(canVar));
      fprintf('\n\nResidual axes  (total = %-3.4f):\n',sum(resVar));
      fprintf('  %3.4f',resVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %3.4f',cumsum(resVar));
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
%
result.SSr        = SSr;
if (partial>0)
   result.SSr_W   = SSr_W;
end
result.SSe        = SSe;
result.SSt        = SSt;

result.MSr        = MSr;
result.MSe        = MSe;

result.X          = X;
result.Y          = NaN;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evects,evals] = sub_decompose(M)
% - eigenvalue decomposition

[evects,evals] = f_eig(M);
evals          = evals'; % row vector of eigenvalues (1 for each axis)

% Only need n-1 axes for n objects:
% evects = evects(:,1:(end-1));
% evals  = evals(1:(end-1));

% Discard eigenvalues = 0:
% tol    = sqrt(eps);
% idx    = find (abs(evals)>tol);
% evects = evects(:,idx);
% evals  = evals(idx);
