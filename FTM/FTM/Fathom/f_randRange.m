function res = f_randRange(minval,maxval,n);
% - returns n random integers ranging from min to max
%
% USAGE: res = f_randRange(min,max,n);
%
% -----Input/Output:-----
% min, max = integers specifying range of random #'s
% n        = number of random #'s to return
% res      = column vector of random integers
%
% SEE ALSO: f_shuffle, f_boot

% -----References:----
% modified after randRange by Ione Fine<ifine@psy.ucsd.edu>, 7/2000 
% http://www-psy.ucsd.edu/~ifine

% -----Author:-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

if (nargin < 3), error('Not enough input parameters');  end;
if (nargin < 4), kind = 1; end; % Use Uniform distribution by default

if (minval>maxval), error('MIN is greater than MAX !'); end;

vec = rand(n,1); % column vector of random #'s

% convert range minval:maxval
vec = ceil((1+maxval-minval)*vec); 
vec(find(vec+1==1)) = 1; % make any zeros into ones 
res = vec+minval-1; 
