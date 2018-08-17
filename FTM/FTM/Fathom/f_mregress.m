function model = f_mregress(x,y,iter,perm,verb)
% - Multiple Linear Regression via Least Squares Estimation
%
% Usage: model = f_mregress(x,y,iter,perm,verb)
%
% x    = matrix of independent variables         (column-wise);
% y    = column vector of dependent variable
% iter = # of iterations for permutation test    (default = 0)
% perm = permute residuals instead of raw data   (default = 1)
% verb = verbose output of results to display    (default = 1)
%
% model = structure with the following fields:
% .F_stat    = F-statistic
% .F_para_p  = parametric  p-value for 1-tailed test of F_stat
% .F_perm_p  = permutation p-value for 1-tailed test of F_stat
% .t_stat    = t-statistic for partial regression coefficients
% .t_para_p  = parametric  p-value for 2-tailed test of t_stat
% .t_perm_p  = permutation p-value for 2-tailed test of t_stat
% .R2        = coefficient of multiple determination (goodness-of-fit)
% .R2adj     = adjusted R2
% .yfit      = fitted values of y
% .b         = regression coefficients (1st value is y-intercept)
% .resid     = residuals

% -----Notes:-----
% This function solves the equation such that:
% y = b(0) + b(1)*(X(:,1)) + b(2)*(X(:,2)) ...+ b(k)*(X(:,k))
% where k = # of predictor variables.
%
% The regression coefficients are computed using Least Squares estimation
% (via the "\" operator), which is preferred over methods that require
% taking the inverse of a matrix. R2, the coefficient of multiple determination,
% is a measure of goodness-of-fit and gives the proportion of variance of
% Y explained by X. The coefficient is adjusted for the number of predictors and
% sample size in R2adj.
%
% Parametric (and optional permutation) tests of significance for the F- and
% t-statistics are performed. The permutation test is conducted when iter > 0
% and allows for permutation of either the raw data or the residuals of the
% full regression model. Permutation of the raw data involves random permutation
% of the rows (= observations) of Y relative to the rows of X. The permutation
% test is preferred over the parametric test when the data are non-normal.
% Permutation of the residuals (vs. the raw data) is preferred when data have
% extreme values (i.e., outliers).
%
% This function has been tested against Legendre & Casgrain's regressn.exe
% program and SAS JMP and gives similar output.

% -----Dependencies:-----
% Calculation of parametric p-values for F and t require fpdf.m and tcdf.m from
% the Matlab Statistics Toolbox, respectively; these could be replaced by df.m
% and dt.m from the free Stixbox Toolbox.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp. (pp.517, 606-612)
% Legendre, P. 2002. Program for multiple linear regression (ordinary or
%   through the origin) with permutation test - User's notes. Depart. of
%   Biological Sciences, University of Montreal. 11 pages. Available from:
%   <http://www.fas.umontreal.ca/biol/legendre/>
% Legendre, P. 2007. Studying beta diversity: ecological variation partitioning
%   by multiple regression and canonical analysis. Journal of Plant Ecology,
%   doi:10.1093/jpe/rtm001
% Neter, J., W. Wasserman, & M. H. Kutner. 1989. Applied linear regression
%   models. 2nd Edition. Richard D. Irwin, Inc. Homewood, IL.

% -----Author:-----
% by David L. Jones, June-2001
% with help from news://comp.soft-sys.matlab
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2002: overhauled permutation test, added calculation of F,t,& R2,
%           & tabular display of results
% Oct-2004: wrap up results into a structure
% Feb-2008: removed nested structures in output; don't transpose t_stat;
%           transpose t_perm_p; changed t-tests from one- to two-tailed
% Mar-2008: added R2adj; preallocated varibles; removed a loop in permutation
%           test; check that Y is a colmn vector

% -----Check input & set defaults:-----
if (nargin < 3), iter = 0; end; % set default iterations to 0
if (nargin < 4), perm = 1; end; % permute residuals by default
if (nargin < 5), verb = 1; end; % set default verbose to 0

% Check that x & y have compatible dimensions:
if (size(x,1) ~= size(y,1)),
   error('X & Y must have the same number of rows!');
end

% Check that y is univariate:
if (size(y,2)>1)
   error('Y must be a column vector specifying a UNIVARIATE response variable!')
end
% -------------------------------------

% Add the y-intercept vector to x (add a column of 1's):
uno  = ones(length(x(:,1)),1);   % unit vector (Neter 6.17)
j    = uno*uno';                 % unit matrix (Neter 6.18)
xOld = x;                        % keep copy of original x for permutation test;
x    = [uno,x];

[n,k] = size(x);  % n = # observations; k = # predictor variables

b     = x\y;      % get regression coefficients via Least Squares Estimation
yfit  = x*b;      % fitted values of y
resid = y - yfit; % residuals

% ========== ANOVA (See Neter et al., 1989:Chapter 7): =========
% Sum-of-Squares Total:               (Neter 7.29)
SSt = y'*y - (1/n)*(y'*j*y);

% Sum-of-Squares Error:               (Neter 7.30)
SSe = y'*y - b'*x'*y;

% Sum-of-Squares of Regression model: (Neter 7.31)
SSr = (b'*x'*y) - (1/n)*(y'*j*y);

% Mean Square Regression model:       (Neter 7.32)
MSr = SSr/(k-1);

% Mean Square Error:                  (Neter 7.33)
MSe = SSe/(n-k);

% F-ratio:                            (for significance tests; Neter 7.34b)
F_stat = MSr/MSe;

% R^2 = Coefficient of Multiple Determination, a goodness-of-fit: (Neter 7.35)
R2 = SSr/SSt; % alternatively, R2 = 1 - SSe/SSt

% R^2adj = Adjusted for # predictors and # samples (L&L,1998 eq.10.20):
R2adj = 1 - ((1-R2)*[(n-1)/(n-k)]);

% Get estimated variance-covariance matrix:        (Neter 7.43)
s2bMatrix = MSe*f_inv((x'*x));
s2b = diag(s2bMatrix); % extract variances for individual coefficients

% t-statistic for partial regression coefficients: (Neter 7.46b)
t_stat = b./sqrt(s2b);

% Parametric p-value for F (one-tailed):
F_para_p = 1 - fcdf(F_stat,k,n-k);

% Parametric p-values for t (two-tailed):
t_para_p = 2 * tcdf(-abs(t_stat), n-1); % after Matlab's ttest.m

% ==============================================================

%-----Permutation Test for F- and t-statistics:-----
rand('state',sum(100*clock)); % set random generator to new state
noCoefs = length(t_stat);     % # coefficients, excluding intercept
if (iter > 0)
   fprintf('\nPermuting the data %d times...\n',iter-1);

   % Preallocate results array:
   randF          = zeros(iter,1);
   randT{noCoefs} = NaN;
   for m=1:noCoefs % variable # of coefficients
      randT{m} = zeros(iter,1);
   end

   % Random permutation:
   for i = 1:(iter-1)
      if (perm>0) % permute residuals:
         residPerm = f_shuffle(resid);
         temp = f_mregress(xOld,residPerm,0,0,0); % need temporary variable since using structures
      else        % permute raw data:
         yPerm = f_shuffle(y,4);               % randomize order of rows only
         temp  = f_mregress(xOld,yPerm,0,0,0); % need temporary variable since using structures
      end
      randF(i) = temp.F_stat;  % keep list of randomized stats (move from structure to vector)
      for m=1:noCoefs
         randT{m}(i) = temp.t_stat(m); % move from  structure to cell array
      end
   end

   % Compute permuted p-values:
   jF = find(randF >= F_stat); % get randomized values >= to F statistic

   jT{noCoefs}        = NaN;
   t_perm_p(noCoefs)  = NaN;
   
   for m = 1:noCoefs
      % Two-tailed test:
      jT{m} = find(randT{m} >= abs(t_stat(m)));  % get randomized values >= to abs(t-stat)
      % Count those vales & convert to two-tailed probability:
      t_perm_p(m) = ((length(jT{m})+1)./iter)*2;
   end

   F_perm_p = (length(jF)+1)./iter; % count those vales & convert to probability
else
   F_perm_p = NaN;
   t_perm_p= repmat(NaN,1,noCoefs);
end
%-------------------------------

if (verb>0)% send output to display:

   fprintf('=====================================================================\n');
   fprintf(' Multiple Linear Regression via QR Factorization:\n');
   fprintf('---------------------------------------------------------------------\n');
   fprintf('R2            R2adj            F-stat        para-p           perm-p \n');
   fprintf('---------------------------------------------------------------------\n');
   fprintf('%-13.5f %-13.5f %-13.5f %-13.5f %-13.5f \n',R2,R2adj,F_stat,F_para_p,F_perm_p);
   fprintf('---------------------------------------------------------------------\n');
   fprintf('\n');
   fprintf('---------------------------------------------------------------------\n');
   fprintf('Variable      b             t-stat        parametric-p  permutation-p\n');
   fprintf('---------------------------------------------------------------------\n');

   for m=1:noCoefs
      if (m==1)
         fprintf('%13s %-13.5f %-13.5f %-13.5f %-13.5f \n',['intercept'],b(m),t_stat(m),t_para_p(m),t_perm_p(m));
      else
         fprintf('%13d %-13.5f %-13.5f %-13.5f %-13.5f \n',(m-1),b(m),t_stat(m),t_para_p(m),t_perm_p(m));
      end
   end

   fprintf('---------------------------------------------------------------------\n\n');
   if (perm>0)
      fprintf('# permutations of residuals = %5d \n',iter-1);
   else
      fprintf('# permutations of rawdata = %5d \n',iter-1);
   end
   fprintf('F-test is one-tailed, t-tests are two-tailed \n');
   fprintf('=====================================================================\n');
end


% Wrap up results into a structure:
model.F_stat   = F_stat;
model.F_para_p = F_para_p;
model.F_perm_p = F_perm_p;
model.t_stat   = t_stat;
model.t_para_p = t_para_p;
model.t_perm_p = t_perm_p';
model.R2       = R2;
model.R2adj    = R2adj;
model.yfit     = yfit;
model.b        = b;
model.resid    = resid;
