function Xr = f_ranging(x,type)
% - scale columns of x to range [0 to 1] or [-1 to 1]
%
% Usage: Xr = f_ranging(x,type);
%
% x    = input matrix
% type = 1:[0 to 1]
%        2:[-1 to 1]
%
% Xr   = scaled matrix
%
% SEE ALSO: f_center, f_stnd, f_transform

% -----Notes:-----
% This function is used to rescale the values of a matrix to range from 0
% to 1 (or -1 to 1), column-wise, following the method of Sneath & Sokal
% (1973). Clarke (1993) suggested that values of species abundance data
% used for SPECIES ordinations or cluster analysis be transformed by the
% total abundance of each species. This can be accomplished with
% f_transform(x,6). F_RANGING(x) also rescales % the data appropriately,
% but may be preferred over F_TRANSFORM as each value can be readily
% interpreted as a proportion of the maximum value for that column (or
% species).

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam. xv + 853 pp.

% -----Author:-----
% by David L. Jones, April-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jul-2005: added code for ranging [-1 to 1]
% Mar-2011: swapped max() and abs() for type==2; added switch-case logic;
%           updated documentation regarding f_transpose

switch type
   case 1
      % eq. 1.11:
      Xr = (x - repmat(min(x),size(x,1),1)) ./ repmat( (max(x)-min(x)) ,size(x,1),1);
   case 2
      % eq. 1.10
      Xr = x ./ repmat(max(abs(x)),size(x,1),1);
   otherwise
      error('TYPE must be 1 or 2!');
end
