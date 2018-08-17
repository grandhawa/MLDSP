function Xstd = f_stnd(x,y)
% - standardize values of a matrix, column-wise (= z-scores)
%
% Usage: Xstd = f_stnd(x,{y});
%
% x    = matrix to standardize
% y    = optional, when present X will be standardized according to Y
%
% Xstd = standardized matrix
%
% SEE ALSO: f_center, f_ranging, f_transform

% -----Notes:-----
% This function is used to standardize the values of a matrix, column-wise,
% creating z-scores. Data matrices containing variables measured in
% different units or on different scales should be standardized before
% being used in ordinations, etc. This makes the variables compatible for
% analyses by putting them on the same scale and giving them equal weight.
% This method is equivalent to Clarke's (1993) 'normalization' of
% environmental variables.
%
% Standardization subtracts the mean and divides by the standard deviation,
% column-wise, so each variable has mean = 0 and stdev = 1;
%
% For ordinations or cluster analyses of VARIABLES (vs. observations) you
% should use F_RANGING instead.

% -----References:-----
% Clarke, K. R. 1993. Non-parametric multivariate analyses of changes
%   in community structure. Aust. J. Ecol. 18: 117-143.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp.

% -----Author:-----
% by David L. Jones, March-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2003: prevent divide by 0
% Jan-2004: stardardize Y based on X

% -----Check input & set defaults:-----
if (nargin<2)
   by_Y = 0; % default don't stardardize by Y
else
   by_Y = 1;
   if size(x,2) ~= size(y,2);
      % error('X & Y must have same # of columns!')
%       OCT-2009: commented out, is this necessary??? [DLJ]
   end
end 

if (by_Y<1)
   
   % -----Center:-----
   Xctr = x - repmat(mean(x),size(x,1),1);
   
   % -----Divide by stardard deviation:-----
   denom = std(x);
   denom(find(denom == 0)) = eps ; % prevent divide by 0 errors:
   Xstd = Xctr ./ repmat(denom,size(x,1),1);
   
else
   
   Xctr = x - repmat(mean(y),size(x,1),1);   
   denom = std(y);
   denom(find(denom == 0)) = eps ;
   Xstd = Xctr ./ repmat(denom,size(x,1),1);
end
