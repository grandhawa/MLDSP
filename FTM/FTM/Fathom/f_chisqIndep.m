function result = f_chisqIndep(Y,iter)
% - Chi-square test of independence (= homogeneity of proportions)
%
% USAGE: result = f_chisqIndep(Y,iter)
%
% Y    = R x C table of observed counts
% iter = number of iterations for permutation test         (default = 0)

% -----Notes:-----
% The data supplied as Y must be counts and NOT proportions.

% -----References:-----
% http://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test
% http://udel.edu/~mcdonald/statchiind.html

% -----Author:-----
% by David L. Jones, Mar-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 2), iter   = 0; end % default iterations for permutation test

% Check size of input:
[r,c] = size(Y);
if (r<2) || (c<2)
   error('Y must have at least 2 rows and 2 columns!')
end

% Check if data are proportions:
if (sum(Y<0)>=1)
   error('Y must be a table of counts NOT proportions!')
end
% -------------------------------------

% Expected frequencies:
E = (repmat(sum(Y,2),1,size(Y,2)) .* repmat(sum(Y),size(Y,1),1) )/sum(Y(:));

% Chi-square statistic:
stat =  sum(sum((((Y - E).^2)./E),2)); %

% Degrees of Freedom:
df = (r-1)*(c-1);

%-----Permutation tests:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   randStat = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      randResult  = f_chisqIndep(f_shuffle(Y,3),0); % permuted Chi-square statistic
      randStat(i) = randResult.stat;                % collect results
   end
   
   j = find(randStat <= stat); % get randomized stats >= to observed statistic
   p = (length(j)+1)./(iter); % count values & convert to probability
else
   p        = NaN;
end
%-----------------------------



% Wrap results up into a structure:
result.stat     = stat; % Chi-square statistic
result.df       = df;   % degrees of freedom
result.p_perm   = p;    % permutation based significance level
result.p_para   = 1 - chi2cdf(stat,df); % parametric p-value

