function Xctr = f_center(x,grp)
% - center data on column mean
%
% USAGE: Xctr = f_center(x);
%
% Xctr  = centered by subtracting mean, column-wise
% grp   = optional column vector of integers specifying group membersip
%
% SEE ALSO: f_stnd, f_ranging, f_transform

% -----Notes:-----
% This function is used to zero center a matrix, column-wise, and should be
% used to eliminate size differences between variables. This is done by
% subtracting the column mean from each observation.
% 
% Providing a grouping variable allows an observation's group mean to be
% subtracted rather than the global mean.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam. page:38

% -----Author(s):-----
% by David L. Jones, March-2003
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% Jan-2011: added support for grouping vector

% -----Check input:-----
if (nargin > 1)
   if (size(x,1) ~= size(grp,1)), error('Size mismatch b/n X & GRP'); end
end
% ----------------------

if (nargin < 2) % no grouping variable
   Xctr  = x - repmat(mean(x),size(x,1),1);

else            % grouping variable
   uGrp = f_unique(grp); % unique groups, unsorted
   nGrp = numel(uGrp);   % # of groups
   Xctr = nan(size(x));  % preallocate
   
   for i = 1:nGrp
      idx = find(grp == uGrp(i)); % index to rows of this group
      Xctr(idx,:)  = x(idx,:) - repmat(mean(x(idx,:)),size(x(idx,:),1),1);
   end
end
