function result = f_CCorA(Y,X,iter,verb)
% - canonical correlation analysis (CCorA)
% 
% USAGE: result = f_CCorA(Y,X,iter,verb);
% 
% Y    = 1st input matrix                          (rows = obs, col = variables)
% X    = 2nd input matrix                          (rows = obs, col = variables)
% iter = # iterations for permutation test                         (default = 0)
% verb = optionally send results to display                        (default = 1)
%
% result = structure of outputs with the following fields:
%  .V        = eigenvectors for Y
%  .U        = eigenvectors for X
%  .evals    = eigenvalues ( = squared canonical correlations, R2)
%  .Cy       = canonical coefficients for Y
%  .Cx       = canonical coefficients for X
%  .yScores  = canonical ordination scores for Y
%  .xScores  = canonical ordination scores for X
%  .xCor     = correlation of canonical axes with X variables
%  .yCor     = correlation of canonical axes with Y variables
%  .trc      = trace statistic
%  .grs      = greatest root statistic
%  .p        = structure of randomized probabilities: p.trace (trace) and
%              p.grs (greatest root)
%  .YX_R2    = R^2 for Y|X
%  .YX_R2adj = adjusted R^2 for Y|X
%  .XY_R2    = R^2 for X|Y
%  .YX_R2adj = adjusted R^2 for X|Y
% 
% SEE ALSO: f_CCorAplot, f_rda, f_cap, f_ncap, f_cda

% -----Notes:-----
% Canonical correlation analysis performs a symmetric analysis using 2
% input matrices, whereas asymmetric analyses (like RDA and multiple
% regression) considers one set of variables is dependent on the other.
% 
% This implementation follows Borcard et al. (2011) and Legendre & Legendre
% (2012) by including the redundancy statistics, which allows one to
% assess whether the correlated variation in the two data sets is large
% with respect to the overall variation.
% 
% CANONICAL COEFFICIENTS: provide the contributions of Y and X to the
% canonical axes.
% 
% Note Legendre & Legendre (unpublished) and vegan's CCorA.R specify using SVD
% vs. eigenanalysis to calculate the eigevectors and eigenvalues. However, this
% function performs an eigenanalysis of the symmetric matrix K'*K (or
% K*K'). From news://comp.soft-sys.matlab: For symmetric, positive definite
% matrices, the eigenvalue and singular value decompositions are
% essentially the same, but the eigenvalue calculation is faster. We use
% SVD in PINV because we are not assuming the matrix is symmetric, positive
% definite, or even square. -- Cleve Moler <moler@mathworks.com>.

% -----References:-----
% Borcard, D., F. Gillet, and P. Legendre. 2011. Numerical Ecology with R.
%  Springer, NY.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% Legendre, P. & L. Legendre. 2012. Numerical ecology. 3rd English ed.
% 
% Equations are from the (unpublished) 3rd edition of Legendre & Legendre (2012).

% -----Author:-----
% by David L. Jones, Mar-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2013: now compatible with changes made to f_pca (Thanks to B. Barnes for
%           identifying this)

% -----Check input & set defaults:-----
if (nargin < 3), iter = 0; end % no permutation test by default
if (nargin < 4), verb = 1; end % send output to display by default

if size(Y,1) ~= size(X,1),
   error('Y & X must have same # rows!')
end
% -------------------------------------

[n,ncY] = size(Y);   % # observations & variables in Y
ncX     = size(X,2); % # variables in X

% Center variables:
Y = f_center(Y); 
X = f_center(X); 

% Replace data by PCA object scores to avoid issues with sparse matrices,
% multicollinearity, etc. (after P. Legendre's CCorA documentation in vegan):
Yold = Y; % keep original data
Xold = X;
Y    = f_pca(Y,0,1,sqrt(eps));
X    = f_pca(X,0,1,sqrt(eps));

% Extract scores:
Y = Y.scores;
X = X.scores;

% Covariance matrices (eq. 11.46):
Syy = (Y'*Y);
Sxx = (X'*X);
Syx = (Y'*X);   % interaction effect
% Sxy = (X'*Y); % interaction effect

% Summarize the correlation structure b/n Y & X:
K = f_inv(chol(Syy')) * Syx * f_inv(chol(Sxx));

% Eigenanalysis:
[V,evals] = f_eig(K*K'); % for Y
U         = f_eig(K'*K); % for X

% Trim null eigenvalues/axes:
s     = min([ncY ncX (n-1)]); % # non-zero canonical eigenvalues
V     = V(:,1:s);
U     = U(:,1:s);
evals = evals(1:s);

% Canonical coefficients:
Cy = f_inv(chol(Syy))*V; % (eq. 11.53)
Cx = f_inv(chol(Sxx))*U; % (eq. 11.54)

% Project data into canonical space:
yScores = Y*Cy; % (eq. 11.55)
xScores = X*Cx; % (eq. 11.56)

% Correlation of original Y variables with each canonical axis (Yold in
% space of Y):
yCor = zeros(ncY,s); % preallocate
for i = 1:s
   for j = 1:ncY
      yCor(j,i) = f_corr(Yold(:,j),yScores(:,i)); % (eq. 11.57)
   end
end

% Correlation of original X variables with each canonical axis (Xold in
% space of X):
xCor = zeros(ncX,s); % preallocate
for i = 1:s
   for j = 1:ncX
      xCor(j,i) = f_corr(Xold(:,j),xScores(:,i)); % (eq. 11.58)
   end
end

% Stats:
trc = sum(evals); % canonical trace statistic
grs = evals(1);   % greatest root statistic

% Redundancy Statistics:
rdaYX    = f_rda(Y,X,0,0,0);
YX_R2    = rdaYX.R2;
YX_R2adj = rdaYX.R2adj;

rdaXY    = f_rda(X,Y,0,0,0);
XY_R2    = rdaXY.R2;
XY_R2adj = rdaXY.R2adj;

%-----Permutation tests:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   randTrc = zeros(iter-1,1); % preallocate
   randGrs = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      
      % Permute obs (rows)
      Yperm = f_shuffle(Y,4);
      
      % Eigenanalysis:
      Kperm            = f_inv(chol((Yperm'*Yperm)'))*(Yperm'*X)*f_inv(chol(Sxx));
      [null,evalsPerm] = f_eig(Kperm*Kperm');

      % Collect statistics:
      randTrc(i) = sum(evalsPerm(1:s)); % permuted trace stat
      randGrs(i) = evalsPerm(1);        % permuted greatest root stat
      
   end
   j1 = find(randTrc >= trc); % get randomized stats >= to observed statistic
   j2 = find(randGrs >= grs);
   
   p.trace = (length(j1)+1)./(iter); % count values & convert to probability
   p.grs   = (length(j2)+1)./(iter);
else
   p.trace = NaN;
   p.grs   = NaN;
end
%-----------------------------


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('      Canonical Correlation Analysis (CCorA:\n');
   fprintf('--------------------------------------------------\n');
   fprintf('Trace Stat    = %-3.4f  p =  %3.5f \n',trc,p.trace);
   fprintf('Greatest Root = %-3.4f  p =  %3.5f \n',grs,p.grs);
   fprintf('No. of permutations = %d \n',iter);
   fprintf('--------------------------------------------------\n');
   fprintf('Canonical Correlations:\n');
   fprintf('  %-3.4f',evals.^0.5);
   fprintf('\n')
   fprintf('Squared Canonical Correlations (= delta^2):\n');
   fprintf('  %-3.4f',evals);
   fprintf('\n--------------------------------------------------\n');
   fprintf('Redundancy Statistics:\n');
   fprintf('Y|X: R2 = %-3.4f   R2adj =  %3.5f \n',YX_R2,YX_R2adj);
   fprintf('X|Y: R2 = %-3.4f   R2adj =  %3.5f \n',XY_R2,XY_R2adj);
   fprintf('==================================================\n');
end


% Wrap results up into a structure:
result.V        = V;       % eigenvectors (= canonical coefficients) for Y
result.U        = U;       % eigenvectors (= canonical coefficients) for X
result.evals    = evals';  % eigenvalues (= squared canonical correlations, R2)
result.Cy       = Cy;      % canonical coefficients for Y
result.Cx       = Cx;      % canonical coefficients for X

result.yScores  = yScores; % canonical ordination scores for Y
result.xScores  = xScores; % canonical ordination scores for X

result.yCor     = yCor;    % correlation of canonical axes with Y variables
result.xCor     = xCor;    % correlation of canonical axes with X variables

result.trc      = trc;     % trace statistic
result.grs      = grs;     % greatest root statistic
result.p        = p;       % structure of randomized probabilities: p.trace (trace)
%                           and p.grs (greatest root)
result.YX_R2    = YX_R2;
result.YX_R2adj = YX_R2adj;
result.XY_R2    = XY_R2;
result.XY_R2adj = XY_R2adj;
