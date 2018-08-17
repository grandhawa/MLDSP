function y = f_issymSim(x,dia)
% - determine if input is square symmetric similarity matrix
%
% USAGE: y = f_isssymSim(x,dia);
%
% x   = input matrix
% dia = test if diagonals are all one (default = 1)
%
% y = 1 if true, 0 if false
%
% SEE ALSO: f_issymdis

% -----Notes:-----
% This function is used to determine if the input matrix X is a square
% symmetric similarity (proximity) matrix having (1) equal number of rows and
% columns, (2) at least 2 rows/columns, (3) identical upper and lower
% tridiagonals, and (4) all 1's along the main diagonal

% -----Author(s):-----
% by David L. Jones, Sep-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), dia = 1; end % default check if main diagonals are all 1

[nr,nc] = size(x);


if (nr~=nc) || (nr<2) % symmetrical size, no row/col vectors
   y = 0;
   return
   
elseif (sum(f_unwrap(x,0) - f_unwrap(x',0))~=0) % symmetrical tridiagonals
   y = 0;
   return
   
elseif (dia>0) && (sum(diag(x)~=1)>0)           % only 1's along diagonal
   y = 0;
   return
   
else
   y = 1;
end

