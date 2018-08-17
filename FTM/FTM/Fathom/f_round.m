function y = f_round(x,n)
% - round data to specified number of decimal places
%
% USAGE: y = f_round(x,n);
%
% x = input vector or matrix
% n = number of decimal places
% y = x rounded to n decimal places

% -----Author:-----
% by David L. Jones, Oct-2005
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

y = round(x*[10^n])/[10^n];
