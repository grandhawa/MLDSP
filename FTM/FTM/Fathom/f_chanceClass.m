function result = f_chanceClass(y,suc,iter,verb)
% - classification success expected by chance (= proportional chance criterion)
%
% USAGE: result = f_chanceClass(y,suc,iter,verb)
%
% y    = column vector of integers specifying group membership
% suc  = observed classification success rate (proportion)        (default = [])
% iter = # iterations for permutation test                         (default = 0)
% verb = optionally send results to display                        (default = 1)
%
% result = structure of results with the following fields:
%   .grp = group
%   .c   = chance probability (by group)
%   .tot = overall chance probability
%   ,p   = randomized probablilty
%
% SEE ALSO: f_errRate

% -----Notes:-----
% The proportional chance criterion is used to determine if the classification
% success rate of a discriminant analysis is better than a null model that
% assigns group membership via random allocation. It provides the probability
% that observations will be classified to the correct group merely by chance and
% is based on the relative proportions of group size.
% 
% Note that when groups sizes are balanced, the overall classification
% success rate expected by chance is 1/g, where g = the number of groups.
% For unbalanced data, the chance classification success rate is greater.
% For highly unbalanced data where one group tends to domminate, the chance
% classification rate increases towards 100%.
% 
% The value 'p' provides the randomized probability that the observed
% overall classification success rate is no better than that expected by
% chance alone.

% -----References:-----
% McGarigal, K., S. Cushman, and S. Strafford. 2000. Multivariate
%   statistics for wildlife and ecology research. Springer-Verlag, New York.
% Morrison, D. G. 1969. On the interpretation of discriminant analysis.
%   Journal of Marketing Research 6:156-163.

% -----Author:-----
% by David L. Jones, Feb-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% May-2012: updated documentation; now supports significance test
% Feb-2013: updated documentation

% -----Set defaults & check input:-----
if (nargin < 2), suc  = []; end % default no success rate provided
if (nargin < 3), iter = 0;  end % no permutation test by default
if (nargin < 4), verb = 1;  end % default send output to display

% Check success rate:
if ~isempty(suc) && (suc<0 || suc>1)
   error('SUC must be between 0 and 1!');
end
% -------------------------------------

% Process groups:
y     = y(:);         % force column vector
n     = size(y,1);    % get # rows
uGrp  = f_unique(y);  % unique groups, unsorted
nGrp  = length(uGrp); % # unique groups
theta = nan(nGrp,1);  % preallocate

% Count the number of observations in each group:
for i=1:nGrp
   theta(i) = sum(y==uGrp(i));
end

% Convert counts to proportions:
theta = theta./sum(theta);

% Proportional chance criterion:
c    = theta.^2;     % by group
tot = sum(theta.^2); % overall

%-----Permutation tests:-----
if (~isempty(suc) && iter>0)
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   rSucc = zeros(iter-1,1); % preallocate
   for i = 1:(iter-1) % observed value is considered a permutation
      rSucc(i) = (sum((y - f_shuffle(y,1)==0))/n); % randomized success rate
   end
   j = find(rSucc >= suc);    % get randomized stats >= to observed statistic
   p = (length(j)+1)./(iter); % count values & convert to probability
else
   p = NaN;
end
%-----------------------------

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('           PROPORTIONAL CHANCE CRITERION:\n'                 );
   fprintf('--------------------------------------------------\n' );
   
   fprintf('Group        Correct  \n');
   for j=1:nGrp
      fprintf('%s %d %s %10.2f %s \n','  ',uGrp(j),'     ',c(j)*100,'%');
   end
   fprintf('\n');
   fprintf('Total Correct = %4.2f %s \n',tot*100,'%');
   fprintf('--------------------------------------------------\n' );
   
   if (~isempty(suc) && iter>0)
      fprintf('\n');
      fprintf('Mean randomized classification success = %4.2f %s \n',mean(rSucc)*100,'%');
      fprintf('p =  %3.5f \n',p');
      fprintf('No. of permutations = %d \n',iter);
      fprintf('--------------------------------------------------\n' );
   end
end
% ---------------------------------

% Wrap results up into a structure:
result.grp = uGrp; % group designation
result.c   = c;    % probability (by group)
result.tot = tot;  % overall probability
result.p   = p;    % randomized significance level
