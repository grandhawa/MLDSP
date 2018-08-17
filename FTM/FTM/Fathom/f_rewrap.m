function xDis = f_rewrap(x)
% - wraps vector into symmetric distance matrix (reverses f_unwrap)
%
% USAGE: xDis = f_rewrap(x,method)
%
% x    = lower tridiagonal extracted from a symmetric distance matrix by f_unwrap
%
% xDis = square symmetric distance matrix
%
% SEE ALSO: f_unwrap, squareform

% -----References:-----
% Legendre, P. & L. Legendre. 2012. Numerical ecology. 3rd English ed.
%   Elsevier Science BV, Amsterdam. [p. 552]
% 
% https://www.mathway.com/

%-----Notes:-----
% This function reverses the effect of f_unwrap by taking a column vector
% defining the lower tridiagonal of a square symmetric distance matrix (usually
% obtained from f_unwrap) and wraps it up into square symmetric form.
%
% This code uses only ONE loop and the algorithm exploits the pattern between
% the size of the lower tridiagonal of a square symmetric matrix (sans the main
% diagonal) and the expected dimensions of the full matrix:
%
% dim of full matrix  = 2 3 4  5  6  7...etc
%                       --------------------
% size of tridiagonal = 1 3 6 10 15 21...etc

% -----Author:-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2003: corrected error in creating upper tri-diagonal
% Apr-2008: use logical indexing vs. find
% Jan-2014: replaced FOR loop with a proper equation to calculate dim

x     = x(:);      % make sure it's a column vector
sizeX = length(x); % size of tridiagonal

% Get dimensions of symmetric matrix:
% There are [n(n-1)/2] values in the lower tridiagonal, so solve for n:
dim = 0.5*(1+sqrt(8*sizeX+1));

% % -----Get dimensions of symmetric matrix:-----
% if (sizeX == 1)
%    dim = 2;
% else
%    dimVar  = 2; sizeVar = 1; % initialize variables
%    while (sizeVar < sizeX);
%       dimVar  = dimVar + 1;  % step thru table shown above
%       sizeVar = (dimVar-1) + sizeVar;
%       if (sizeVar == sizeX)
%          dim = dimVar;
%       elseif(sizeVar > sizeX)
%          error('Input vector wrong size !');
%       end
%    end
% end
% % ---------------------------------------------

xDis     = zeros(dim,dim);              % preallocate results matrix
ai       = logical(tril(ones(dim),-1)); % get indices for lower diagonal
xDis(ai) = x;                           % fill lower tridiagonal
xDis     = xDis + xDis';                % fill upper tridiagonal


