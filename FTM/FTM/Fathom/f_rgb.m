function res = f_rgb(n)
% - utility program for selecting color of plot symbols
%
% USAGE res = f_rgb(n)
%
% n   = interger value selecting color
% res = symbol specifying rgb triplet for plots
%
% SEE ALSO: f_symb, f_style

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2003: when used in conjunction with f_symbol, the following was added to
%           prevent duplicate symbol colors when the # of combinations < 101.
% Feb-2008: changed & to &&

% Aug-2009: changed order of first 3 colors


if (n>10 && n<21),  n = n + 1; end;
if (n>20 && n<31),  n = n + 2; end;
if (n>30 && n<41),  n = n + 3; end;
if (n>40 && n<51),  n = n + 4; end;
if (n>50 && n<61),  n = n + 5; end;
if (n>60 && n<71),  n = n + 6; end;
if (n>70 && n<81),  n = n + 7; end;
if (n>80 && n<91),  n = n + 8; end;
if (n>90 && n<101), n = n + 9; end;


while (n>10)
   n = n - 10;
end

% Build cell-array of colors to choose from:
colors{1}  = [1 0 0];          % red
colors{2}  = [0 0.5 0];        % dark green
colors{3}  = [0 0 1];          % blue
colors{4}  = [0.28 0.73 0.94]; % sky
colors{5}  = [0 0 0];          % black
colors{6}  = [0 0 0.5];        % navy
colors{7}  = [0 1 0];          % green
colors{8}  = [1 0.5 0];        % orange
colors{9}  = [0.5 0 0.5];      % purple
colors{10} = [1 0 1];          % magenta

res = colors{n};
