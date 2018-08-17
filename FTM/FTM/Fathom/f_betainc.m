function Y = f_betainc(X,Z,W,tail)
% - incomplete beta function ('non-regularized')
% 
% USAGE Y = f_betainc(X,Z,W,tail);
% 
% X    = defines interval of the integral
% Z,W  = parameters of incomplete beta function
% tail = as 'lower' (= default) or 'upper' 
% 
% SEE ALSO: betainc

% -----Notes:-----
% The function computes the incomplete beta function and differes from the
% Matlab BETAINC function which computes the 'regularized' incomplete
% beta function (i.e., the cumulative beta distribution). The regularized
% incomplete beta function multiplied by the complete beta function yields the
% actual 'incomplete beta function'.

% -----References:-----
% http://en.wikipedia.org/wiki/Beta_function:

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 4), tail = 'lower'; end % default lower tail
% -------------------------------------

% Get incomplete beta function:
Y = betainc(X,Z,W,tail) .* beta(Z,W);
