function result = f_CCorA_PL(Y,X,iter,verb)
% - CCorA based on P. Legendre's vegan code
%
% USAGE: result = f_CCorA_PL(Y,X,iter,verb);
%
% Y    = 1st input matrix                          (rows = obs, col = variables)
% X    = 2nd input matrix                          (rows = obs, col = variables)
% iter = # iterations for permutation test                         (default = 0)
% verb = optionally send results to display                        (default = 1)

%
% result = structure of outputs with the following fields:
%  .U       = eigenvectors (= canonical coefficients) for Y
%  .V       = eigenvectors (= canonical coefficients) for X
%  .evals   = eigenvalues (= squared canonical correlation, R2)
%  .yScores = canonical ordination scores for Y
%  .xScores = canonical ordination scores for X
%  .yCor    =  correlation of canonical axes with Y variables
%  .xCor    =  correlation of canonical axes with X variables
%  .trc     = trace statistic
%  .grs     = greatest root statistic
%  .p       = structure of randomized probabilities: p.trace (trace) and
%             p.grs (greatest root)
%
% SEE ALSO: f_CCorAplot, f_cap, f_ncap, f_rda, f_cda

% -----Notes:-----
% Canonical correlation analysis performs a symmetric analysis using 2
% input matrices, whereas asymmetric analyses (like RDA and multiple
% regression) considers one set of variables is dependent on the other.
%
% The canonical eigenvalues are the squared canonical correlations, so
% yEvals and xEvals should be equivalent.

% -----References:-----
% Anderson, M. J. & T. J. Willis. 2003. Canonical analysis of principal
%   coordinates: a useful method of constrained ordination for ecology. Ecology
%   84(2): 511-525.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. pages 612-616
% 
% Portions of this function are based on ideas implemented in Legendre's
% 'CCorA.R' code in the vegan package for R.

% -----Author:-----
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 3), iter = 0;  end % no permutation test by default
if (nargin < 4), verb = 1;  end % send output to display by default

if size(Y,1) ~= size(X,1),
   error('Y & X must have same # rows!')
end
% -------------------------------------

[n,ncY] = size(Y);   % # observations & variables in Y
ncX     = size(X,2); % # variables in X

% Center data:
Y = f_center(Y);
X = f_center(X);

% Replace data by PCA object scores to avoid issues with sparse matrices,
% multicollinearity, etc. (P. Legendre in vegan's CCorA documentation):
Yold = Y; % keep original data
Xold = X;
Y    = f_pca(Y);
X    = f_pca(X);

% Covariance matrices (eq. 11.21):
Syy = (Y'*Y)/(n-1);
Sxx = (X'*X)/(n-1);
Syx = (Y'*X)/(n-1); % interaction effect

% Singular Value Decomposition (after CCorA.R):
[U,null,V] = f_svd(f_inv(chol(Syy))' * Syx * f_inv(chol(Sxx)));

% Biplot scores (after CCorA.R):
yScores = Y*(f_inv(chol(Syy))*U);
xScores = X*(f_inv(chol(Sxx))*V);

% Correlation of original Y variables with each canonical axis:
yCor = zeros(ncY,size(yScores,2)); % preallocate
for i = 1:s
   for j = 1:ncY
      yCor(j,i) = f_corr(Yold(:,j),yScores(:,i));
   end
end

% Correlation of original Y variables with each canonical axis:
xCor = zeros(ncX,size(xScores,2)); % preallocate
for i = 1:s
   for j = 1:ncX
      xCor(j,i) = f_corr(Xold(:,j),xScores(:,i));
   end
end

% Stats:
trc = sum(evals); % canonical trace statistic
grs = evals(1);   % greatest root statistic


%-----Permutation tests:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   randTrc = zeros(iter-1,1); % preallocate
   randGrs = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      
      % Permute obs (rows)
      Yperm = f_shuffle(Y,4);
      
      % Eigenanalysis:
      [null,evalsPerm] = f_eig((Yperm'*X)*f_inv(Sxx)*(X'*Yperm)*f_inv((Yperm'*Yperm)));
      
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
   fprintf('\n==================================================\n');
end


% Wrap results up into a structure:
result.U       = U;      % canonical coefficients (= eigenvectors) for Y
result.V       = V;      % canonical coefficients (= eigenvectors) for X
result.evals   = evals'; % canonical correlations (= eigenvalues = R2)

result.yScores = yScores; % canonical ordination scores for Y
result.xScores = xScores; % canonical ordination scores for X

result.yCor    = yCor;    % correlation of canonical axes with original Y variables
result.xCor    = xCor;    % correlation of canonical axes with original X variables

result.trc     = trc;     % trace statistic
result.grs     = grs;     % greatest root statistic
result.p       = p;       % structure of randomized probabilities: p.trace (trace)
%                           and p.grs (greatest root)


