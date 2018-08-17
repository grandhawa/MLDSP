function result = f_npManovaPW(yDis,x,iter,verb)
% - a posteriori, multiple-comparison tests
%
% USAGE: result = f_npManovaPW(yDis,x,iter,verb);
%
% yDis = square symmetric distance matrix derived from response variables
% x    = matrix of integers specifying factor levels for objects
%        in yDis (column-wise)
% iter = # iterations for permutation test  (default = 0)
% verb = optionally send results to display (default = 1)
%
% result   = structure with the following fields:
%  .t        = t-statistic
%  .pairList = list of the pair-wise treatment level tests
%  .p        = permutation-based significance (unadjusted)
%  .p_bon    = Bonferonni-adjusted p-values
%  .p_ds     = Dunn-Sidak adjusted p-values
%  .p_holm   = Holms adjusted p-values
% 
% See also: f_npManova, f_adjustP

% -----Notes:-----
% This function is used to perform _a posteriori_, multiple-comparison
% (pair-wise) tests after finding a significant factor effect 
% using 'f_npManova'. This essentially involves performing a number of 
% a single classification (M)ANOVAs using all possible pairs of
% treatment levels from the ANOVA factor specified. The t-statistic
% returned is simply the square-root of the usual F-statistic.
%
% Remember you can make additional, unplanned comparisons as well by
% recoding the treatment levels. For example, you found a significant
% treatment effect which had 4 levels (level 1 was the control). You can
% recode the levels via 'x(find(x>1)) = 2' and simply test the control
% vs. the non-control's, etc.
%
% Three methods of correction for multiple comparisons are provided. The
% Bonferroni and Dunn-Sidak are overly conservative, especially considering
% that at alpha = 0.05 only 1 out of 20 tests will be found to be
% significant by random chance (Anderson, 2000). The Holms adjustment is
% the most powerful of the three methods as it results in rejection of the
% null hypothesis more often (Legendre & Legendre, 2012).

% ----- References:-----
% Anderson, M. J. 2000. NPMANOVA: a FORTRAN computer program for non-parametric
%   multivariate analysis of variance (for any two-factor ANOVA design) using
%   permutation tests. Dept. of Statistics, University of Auckland.
%   (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% 
% Legendre, P. & L. Legendre. 2012. Numerical ecology. 3rd English ed.
%   Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, Nov-2002
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2007: preallocated t,p; wrapped up results into a structure
% Mar-2010: show which test + # iterations in command window, change output when
%           verb=1
% Oct-2012: p-values now adjusted using f_adjustP; all p-values are provided in
%           verbose output
% Mar-2013: updated documentation

% ----Check input & set defaults:-----
if (nargin < 3), iter  = 0; end % default iterations for permutation test
if (nargin < 4), verb  = 1; end % don't send output to display by default

% Extract all pair-wise subsets:
[sDis,sX,pairList] = f_subsetDisPW(yDis,x);

% Get # of pairwise tests:
noTests = size(pairList,1);

% 1-way ANOVA for each pair of treatment levels:
t = repmat(NaN,noTests,1); % preallocate
p = repmat(NaN,noTests,1);
for i=1:noTests
   if (verb>0)
      fprintf('\nPermuting the data %d times for test %d of %d...',iter-1,i,noTests);
   end
   temp = f_npManova(sDis{i},sX{i},iter,0);
	t(i) = sqrt(temp(1).F); % since temp(2:3).F = NaN's  [t is sqrt of F]
	p(i) = temp(1).p;       % since temp(2:3).p = NaN's
end

% Adjust p-values for multiple-comparisons:
p_bon  = f_adjustP(p,'bon');
p_ds   = f_adjustP(p,'ds');
p_holm = f_adjustP(p,'holm');

% Send output to display:
if (verb>0)
   fprintf('\n\n')
   fprintf( '----------------------------------------------------------\n')
   fprintf('Results of pair-wise comparisons between each factor level:\n')
   fprintf( '==========================================================\n')
   fprintf('         t:     p:     p_bon: p_ds:  p_holm:\n')
   for i = 1:size(pairList,1)
      fprintf('%d vs. %d: %6.4f %6.4f %6.4f %6.4f %6.4f \n',pairList(i,1),pairList(i,2),...
         t(i),p(i),p_bon(i),p_ds(i),p_holm(i));
   end
   fprintf( '----------------------------------------------------------\n')
   fprintf('t      = t-statistic \n');
   fprintf('p      = unadjusted p-value \n');
   fprintf('p_bon  = Bonferroni adjusted p-value \n');
   fprintf('p_ds   = Dunn-Sidak adjusted p-value \n');
   fprintf('p_holm = Holms adjusted p-value \n');
   fprintf('\n\n')
end

% Wrap up results into a structure: 
result.t        = t;
result.pairList = pairList;
result.p        = p;
result.p_bon    = p_bon;
result.p_ds     = p_ds;
result.p_holm   = p_holm;
