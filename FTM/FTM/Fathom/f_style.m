function res = f_style(n)
% - utility program for selecting line styles
%
% USAGE res = f_style(n)
%
% n   = interger value selecting
% res = symbol specifying linespec for plots
%
% SEE ALSO: f_rgb, f_symb

% -----Author:-----
% by David L. Jones, Nov-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

while (n>4)
   n = n - 4;
end

% Styles to choose from:
styles = {'-','--',':','-.'};

res = styles{n};