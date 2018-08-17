function s = f_corrSign(x,y)
% - determine if a correlation b/n 2 vectors is positive or negative
%
% USAGE: [r,p] = f_corrSign(x,y);
%
% x,y  = input vectors
% s    = positive (= 1) or negative (-1) correlation
%
% SEE ALSO: f_corr

% -----Author(s):-----
% by David L. Jones, Jan-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
% make sure they're column vectors
x = x(:);
y = y(:);

n = length(x);

if (n ~= length(y))
   error('Input vectors must same size!');
end
% -------------------------------------

% Calculate regression:
model = f_mregress(x,y,0,0,0);

% Determine the sign of the slope:
if model.b(2) < 0
   s = -1;
else
   s = 1;
end
