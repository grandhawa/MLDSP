function cv = f_cv(x)
% - coefficient of variation
% 
% USAGE: result = f_cv(x);
% 
% x  = input matrix (row = observations, col = variables)
% cv = coefficient-of-variation of corresponding cols

% -----Notes:-----
% The coefficient-of-variation is a normalized measured of dispersion defined as
% the ratio of the standard deviation to the mean. It is usually presented as
% the ratio multiplied by 100. It describes the dispersion of a single variable
% independent of it measurement units. Higher values correspond to greater
% dispersion in the variable.
% 
% It is only defined for data with a non-zero mean and most useful for variables
% that are always positive.

% -----Author:-----
% by David L. Jones,<djones@rsmas.miami.edu> Sep-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Get # rows of x:
r = size(x,1);

% Make sure data are all > 0
x = x - repmat(min(x),r,1);

% Calculate coefficient:
cv = (std(x)./mean(x)) * 100;
