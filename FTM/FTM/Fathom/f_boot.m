function [B,T] = f_boot(X,c,n)
% - bootstrap resampling with replacement
%
% USAGE: [B,T] = f_boot(X,c,n);
%
% X = input matrix                    (rows = obs, cols = variables)
% c = resample each column separately (default = 0);
% n = number of bootstrapped samples  (default = same as input X)
%
% B = bootstrapped sample
% T = elements of X not in B
%
% SEE ALSO: f_bootCI, f_grpBoot, f_shuffle, f_randRange

% -----Author:-----
% by David L. Jones, Dec-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2004: now uses unidrnd vs. f_randRange, added output of T
% Jan-2010: replaced '&' with '&&'
% Oct-2010: now works with matrices, can specify n, optional resample columns
%           separately
% Feb-2014: replaced repmat with nan; generate T only when needed

% -----Set defaults and check input:-----
if (size(X,1)==1), X = X(:);    end % handle row vectors
if (nargin < 2), c = 0;         end % default don't resample cols separately
if (nargin < 3), n = size(X,1); end % default same # of samples as X

% Generate T only when necessary:
if (nargout <2), getT = 0; else getT=1; end
% -----------------------------------------

[nr,nc] = size(X);
B       = nan(n,nc); % preallocate

if (c>0) % Resample each column separtely:
   for i = 1:nc
      B(:,i) = X(unidrnd(nr,n,1),i);
   end
   T = NaN; % not applicable here
else % Resample obs with replacement:
   B = X(unidrnd(nr,n,1),:);
   if (getT==1)
      T = setdiff(X,B,'rows'); % elements not resampled
   else
      T = NaN;
   end
end
