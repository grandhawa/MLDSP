function result = f_bootCI(X,conf,iter)
% - bootstrapped confidence interval of the mean for univariate data
%
% USAGE: result = f_bootCI(X,conf,iter);
%
% X     = input data (rows are obs, columns are variables
% conf  = confidence interval                                   (default = 0.95)
% iter  = # iterations for bootstrap resampling                    (default = 0)
% 
% result = structure of results with the following fields:
%  .mean = mean values 
%  .cv   = critical values
%  .conf = confidence level
% 
% SEE ALSO: f_boot

% -----Notes:-----
% For bootstrapped confidence intervals it is recommened that ITER be at least
% 1000.

% -----Author:-----
% by David L. Jones, Sep-2013
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 2), conf  = 0.95; end; % default 95% confidence interval
if (nargin < 3), iter  = 0;    end; % default 0 bootstrap samples

if (conf <= 0) || (conf >= 1)
   error('CI must be >0 and <1)');
end
% ---------------------------------------

nc        = size(X,2);    % get # rows, cols of input x
boot      = NaN(iter,nc); % preallocate result array
boot(1,:) = mean(X);      % get observed mean for each variable
   
% Get bootstrapped means:
for i=2:iter
   boot(i,:) = mean(f_boot(X,1)); % get mean of bootstrapped sample
end

% Bootstrapped critical value:
cv = prctile(boot,[1-(1 - conf)]*100);

diffVar = abs(cv - boot(1,:));


% Wrap results up into a structure:
result.mean = boot(1,:);
result.ci_u = boot(1,:) + diffVar;
result.ci_l = boot(1,:) - diffVar;
result.conf = conf;
