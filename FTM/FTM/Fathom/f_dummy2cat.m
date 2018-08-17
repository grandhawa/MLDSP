function  Y = f_dummy2cat(X)
% - convert dummy codes to a single categorical variable
%
% USAGE: Y = f_dummy2cat(X)
%
% X = binary dummy codes created by f_dummy
% Y = column vector of integers specifying factor levels or group membership
% 
% SEE ALSO: f_dummy

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


% -----Notes:-----
% The purpose of this code is to essentially undo the effects of f_dummy

[nr,nc] = size(X);
Y       = zeros(nr,1); % preallocate

for i = 1:nc
   Y(X(:,i)==1) = i;
end








