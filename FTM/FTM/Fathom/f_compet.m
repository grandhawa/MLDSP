function c = f_compet(x)
% - get index of column containing the maximum value of each row
%
% USAGE: c = f_compet(x);
%
% x = input matrix (row = obserations, columns = classes)
% c = column having maximum value for each row; in the case of tie the
%     first column having the maximum value will be returned

% -----Notes:-----
% This function was written primarily to determine the 'winner' from the
% output of a neural network.

% -----Author:-----
% by David L. Jones, Feb-2004
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

[null,c] = max(x');
c        = c';
