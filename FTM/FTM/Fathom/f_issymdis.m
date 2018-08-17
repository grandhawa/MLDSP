function y = f_issymdis(x,dia)
% - determine if input is square symmetric distance matrix
%
% USAGE: y = f_isssymdis(x,dia);
%
% x   = input matrix
% dia = test if diagonals are all zero (default = 1)
%
% y = 1 if true, 0 if false
% 
% SEE ALSO: f_issymSim

% -----Notes:-----
% This function is used to determine if the input matrix X is a square
% symmetric distance matrix having (1) equal number of rows and columns, (2) at
% least 2 rows/columns, (3) identical upper and lower tridiagonals, and (4) all
% 0's along the main diagonal

% -----Author(s):-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2003: more efficient coding: doesn't take diag(X) unless X is square,
%           looks for nonzero values along diagonal vs. sum(diag(x)).
% Apr-2008: fixed spelling; replaced | with ||; slight rewrite; testing for only
%           zeros along the main diagonal is now optional

% -----Check input & set defaults:-----
if (nargin < 2), dia = 1; end % default check if main diagonals are all 0

[nr,nc] = size(x);


if (nr~=nc) || (nr<2) % symmetrical size, no row/col vectors
   y = 0;
   return

elseif (sum(f_unwrap(x,0) - f_unwrap(x',0))~=0) % symmetrical tridiagonals
   y = 0;
   return

elseif (dia>0) && (sum(diag(x)~=0)>0)           % only 0's along diagonal
   y = 0;
   return

else
   y = 1;
end
