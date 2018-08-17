function Xp = f_perturb(X)
% - randomly perturb the values of a matrix, separately for each column
% 
% USAGE: Xp = f_perturb(X);
% 
% X  = input matrix (rows = observations, cols = variables)
% Xp = perturbed values of X

% -----Notes:-----
% This method perturbs the input matrix X with values having the same properties
% as the original data.

% -----References:-----
% Silverman, B. W. 1986. Density estimation for statistics and data analysis.
%  Chapman and Hall, London (Section 6.4, page 141)

% -----Author:-----
% by David L. Jones, Jun-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

[n,nc] = size(X);   % get size of input
Xp     = nan(n,nc); % preallocate

% Perturb each column separately:
for i=1:nc
   % eq. 3.28 in Silverman (1986):
   h = 1.06 * std(X(:,i)) * n^(-1/5);
   
   % eq. 6.16 in Silverman (1986):
   x_bar   = mean(X(:,i));
   var_x   = var(X(:,i));
   var_k   = 1;
   Xp(:,i) = x_bar + (X(:,i) - x_bar + randn(n,1)*h)/sqrt(1 + h^2*var_k/var_x);
end
