function res = f_symb(n)
% - utility program for selecting plot symbols
%
% USAGE res = f_symb(n)
%
% n   = interger value selecting
% res = symbol specifying linespec for plots
%
% SEE ALSO: f_rgb, f_style

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

while (n>10)
   n = n - 10;
end

% symbols to choose from:
symbols = {'o','+','^','*','p','s','d','.','x','h'};

res = symbols{n};