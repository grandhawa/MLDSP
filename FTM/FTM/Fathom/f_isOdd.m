function sol = f_isOdd(x)
% - determine if integer is odd or even
%
% USAGE: sol = f_isOdd(x);
%
% x   = scalar, vector, or matrix
% sol = boolean (logical) scalar, vector, or matrix
%       indicating 1 = even, 0 = odd

% -----Notes:-----
% Here's an example of how to extract the odd elements from a vector:
% test = [-5:5]
% >> test(f_isOdd(test))
% ans = -5    -3    -1     1     3     5

% -----Author:-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: updated documentation, added to FATHOM

sol = mod(x,2);
sol = logical(sol == 1);
