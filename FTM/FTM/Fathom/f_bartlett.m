function [pval, chisq, df] = f_bartlett(x,groups)
% - Bartlett's test for Homogeneity of Variances
%
% Usage:  [pval, chisq, df] = f_bartlett(x,groups)
%
% x      = column vector of input values
% groups = column vector specifying group membership
%
% pval  = significance level
% chisq = Chi-square statistic
% df    = degrees of freedom

% Under the null of equal variances, the test statistic CHISQ
% approximately follows a chi-square distribution with DF degrees of
% freedom; PVAL is the p-value (1 minus the CDF of this distribution at
% chisq) of the test.
%

%-----Author:-----
% original Author:  KH <Kurt.Hornik@ci.tuwien.ac.at>
% ported to Matlab by Dave Jones, Dec-2001
% originally bartlett_test.m in Octave
% output verified  with "bartlett.sas" in SAS
%
% Nov-2002: modified by David L. Jones to handle input vector X
%           with groups defined by column vector GROUPS,
%           with additional checking of input
%
% This file is part of the FATHOM Toolbox for Matlab.

% -----Check input:-----
if (size(x,2) > 1)
	error('X must be a column vector')
end

if (size(x,1) ~= size(groups,1))
	error('X and GROUPS must have SAME # or rows')
end 
% ----------------------

subGroups = unique(groups);    % unique groups [DJ]
k         = size(subGroups,1); % no of unique groups [DJ]

f = zeros(k, 1); % preallocate
v = zeros(k, 1);

for i = 1:k;
	subX = x(find(groups == subGroups(i)));
	f(i) = length(subX) - 1;
   v(i) = var(subX);
end

f_tot = sum(f);
v_tot = sum(f .* v) / f_tot;
c     = 1 + (sum(1 ./ f) - 1 / f_tot) / (3 * (k - 1));
chisq = (f_tot * log(v_tot) - sum (f .* log(v))) / c;
df    = k-1; % this was k in Octave version
pval  = 1 - chi2cdf(chisq, df);

