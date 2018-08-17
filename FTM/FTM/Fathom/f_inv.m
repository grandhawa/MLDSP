function invX = f_inv(X)
% - matrix inversion via "\" (left division)
%
% USAGE: invX = f_inv(X);
%
% X    = input matrix, preferably a square symmetric matrix
% invX = generalized inverse of matrix X
%
% SEE ALSO: inv, pinv, qr, slash, \, mldivide

% -----References:-----
% See Matlab online documentation for "\"

% -----Author:-----
% by David L. Jones, Oct-2002
% 
% This file is part of the FATHOM Toolbox for Matlab and 
% is released under the GNU General Public License, version 2.

invX = X\eye(size(X));
