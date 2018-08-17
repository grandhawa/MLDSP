function sub = f_randSub(X,n)
% - extract a random subset of N rows from matrix X
%
% USAGE: sub = f_randSub(X,n);
%
% X   = input matrix (or column vector)
% n   = # of rows of X to extract
%
% sub = random subset of rows of X
%
% SEE ALSO: f_boot, f_randRange, f_shuffle

%-----Notes:-----
% This function is used to extract a random subset of the ROWS of input
% matrix X, while maintaining all COLUMNS. This is random sampling WITHOUT
% replacement. This function may be useful in those circumstances when you
% need to randomly sample among groups to obtain a balanced design in
% MANOVA, CDA, etc. 

% -----Author:-----
% by David L. Jones, Aug-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
noRows = size(X,1);
if (n>noRows)
   error('N is greater than the # of rows in X!')
end
% ----------------------

idx = f_shuffle([1:noRows]'); % shuffle row indices
idx = idx(1:n);               % subset shuffled indices
idx = sort(idx);              % re-sort to original order
sub = X(idx,:);               % subset rows

