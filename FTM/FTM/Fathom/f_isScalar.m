function sol = f_isScalar(x)
% - determine if input is a scalar
%
% USAGE: sol = f_isScalar(x);
%
% x   = input
% sol = boolean (logical) indicating 1 = true, 0 = false
%
% SEE ALSO: f_isOdd


% -----Author:-----
% by David L. Jones, Feb-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: replaced & with &&

if (size(x,1)==1) && (size(x,2)==1)
   sol = 1;
else
   sol = 0;
end


