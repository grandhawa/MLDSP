function result = f_ancovaPW(Y,grp,W,iter,verb,method)
% - a posteriori, multiple-comparison tests for Homogeneous Slopes in ANCOVA
%
% USAGE: result = f_ancovaPW(Y,grps,W,iter,verb,method);
%
% Y      = column vector of response variable
% grps   = column vector of whole numbers specifying group membership     (= Rx)
% W      = column vector specifying covariate
% iter   = # iterations for permutation test                       (default = 0)
% verb   = optionally send results to display                      (default = 1)
% method = adjust p-values via bon' (= Bonferroni), 'ds' (= Dunn-Sidak),
%          'holm' (= Holmes), or 'none' (= default no correction applied)
%
% result   = structure with the following fields:
%  .t        = t-statistic
%  .p        = permutation-based significance
%  .pairList = list of the pair-wise treatment level tests
%  .corrMC   = corrections for multiple-comparison tests 
%             (at alpha = 0.05) via Dunn-Sidak & Bonferonni Methods
%
% See also: f_ancova, f_npManovaPW

% -----Notes:-----
% This function is used to perform _a posteriori_, multiple-comparison
% (pair-wise) tests after finding a significant "Treatment by Covariate"
% Interaction Effect using 'f_ancova' (i.e., heterogeneous slopes of the
% regression lines). This essentially involves performing a number of
% 1-way ANCOVAs using all possible pairs of treatment levels from the
% ANCOVA factor specified.
%
% Remember you can make additional, unplanned comparisons as well by
% recoding the treatment levels. For example, if you found a significant
% interaction effect which had 4 levels (level 1 was the control). You can
% recode the levels via 'x(find(x>1)) = 2' and simply test the control
% vs. the non-control's, etc.
%
% Two methods of correction for multiple comparisons are provided. These
% are highly conservative, especially considering that at alpha = 0.05
% only 1 out of 20 tests will be found to be significant by random chance
% (Anderson, 2000).

% -----Author:-----
% by David L. Jones, Jul-2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% ----Check input & set defaults:-----
if (nargin < 4), iter   = 0; end % no permutation test by default
if (nargin < 5), verb   = 1; end % send output to display by default
if (nargin < 6), method = 'none'; end % send output to display by default

% Internal flag specifying a "PW slope run":
PW = 1;

% Force column vectors:
Y   = Y(:);
grp = grp(:);
W   = W(:);

% Extract all pair-wise subsets for Y:
[sY,sGrp,pairList] = f_subsetPW(Y,grp);

% Extract all pair-wise subsets for W:
sW = f_subsetPW(W,grp);

% Get # of pairwise tests:
noTests = size(pairList,1);

% 1-way ANCOVA for each pair of treatment levels:
F_slope = repmat(NaN,noTests,1); % preallocate
p_slope = repmat(NaN,noTests,1);
for i=1:noTests
   if (verb>0)
      fprintf('\nPermuting the data %d times for test %d of %d...',iter-1,i,noTests);
   end
   temp       = f_ancova(sY{i},sGrp{i},sW{i},iter,0,0,PW);
   F_slope(i) = temp.F_slope;
	p_slope(i) = temp.p_slope;
   clear temp;
end

% Adjust p-values:
p_slope = f_adjustP(p_slope,method);

% Send output to display:
if (verb>0)
   fprintf('\n\nPair-wise comparisons of SLOPES between each treatment level:')
   fprintf(  '\n----------------------------------------------------------')
   for i = 1:size(pairList,1)
      fprintf('\n%d vs. %d: F = %-7.4f p = %-5.4f',pairList(i,1),pairList(i,2),F_slope(i),p_slope(i));
   end
   fprintf('\n\n')
   fprintf('Correction method: %s\n',method);
end

% Wrap up results into a structure: 
result.F_slope  = F_slope;
result.p_slope  = p_slope;
result.pairList = pairList;
