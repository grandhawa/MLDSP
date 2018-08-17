function se = f_stdErr(x)
% - returns the standard error
%
% USAGE: se = f_stdErr(X);
%
% x  = column vector or matrix of input data
% 
% se = standard error of each variable in X (column-wise)
%
% SEE ALSO: var, std

% -----Author:-----
% by David L. Jones, June-2007
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2009: corrected documentation

nr = size(x,1);
se = sqrt(var(x)/nr);
