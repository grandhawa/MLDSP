function Y = f_diag(X,A)
% - replace diagonals of a square matrix
% 
% USAGE: Y = f_diag(X,A)
% 
% X = square input matrix
% A = scalar value to replace diagonals of X matrix with
% 
% Y = new matrix with A's along the diagonal

% -----Author:-----
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
[r,c] = size(X);
if (r~=c), error('X must be a square matrix!'); end

if (numel(A)>1), error('A must be a scalar!'); end
% ----------------------   

% Put value of A along X's diagonal:
X(eye(size(X))==1) = A;

% Rename for output:
Y = X;
