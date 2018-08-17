function minMax = f_range(x)
% - return the min and max values of a vector
%
% USAGE: minMax = f_range(x);
%
% x = row or column vector

% -----Author:-----
% by David L. Jones, Dec-2002
% with help from news://comp.soft-sys.matlab

%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

[nr,nc] = size(x);

if (nr ~= 1) & (nc ~= 1)
	error('X must be a row or column vector');
end

x = x(:);
minMax = [min(x) max(x)];

