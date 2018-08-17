function y = f_unwrap(x,test)
% - unwrap lower tri-diagonal (w/o diag) of symmetric distance matrix into a column vector
%
% USAGE: y = f_unwrap(x);
%
% x = symmetric distance matrix
% y = column vector of extracted elements
%
% SEE ALSO: f_rewrap

% -----References:-----
% after USENET article by Eugene Gall<eugenegall@aol.com>
% posted to news://comp.soft-sys.matlab

% -----Authors:-----
% by David L. Jones
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% note the use of "triu" to return lower tridiagonal (DJ)

% Mar-2002: added checking of input
% Apr-2008: use logical indexing vs. find; added internal flag 'test' for
%           compatibility with f_issymdis

% -----Set defaults & check input:-----
if (nargin < 2), test = 1; end % internal flag, test f_issymdis by default

if test>0
   if (f_issymdis(x) == 0)
      error('Requires square symmetric distance matrix');
   end
end

y = x(~triu(ones(size(x))));
