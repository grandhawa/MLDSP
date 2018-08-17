function txt = f_num2cell(X)
% - convert a vector of numbers to a cell array of strings
% 
% USAGE: txt = f_num2cell(X)
% 
% X   = input vector of numbers
% txt = cell array of strings

% -----Author:-----
% by David L. Jones, Dec-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

txt = cellstr(num2str(X(:))); % create labels from numbers
