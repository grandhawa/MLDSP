function gRel = f_grpRel(x,grps)
% - returns the relative proportion of X (column-wise) separately for groups
%
% USAGE: gRel = f_grpRel(x,grps);
%
% x    = column vector or matrix of input data
% grps = column vector of integers specifying group memberhip
%
% gRel = relative proportion of each group (row-wise) for each variable in X (column-wise)
%
% SEE ALSO: f_grpMean

% -----Notes:-----
% This function is used to calculate the relative proportion of each
% variable in matrix X (column-wise) separately for each group specified by
% GRPS. The number of rows in gRel will equal the number of unique groups
% specified in GRPS. The values in gRel are sorted row-wise in ascending
% order according to the groups specifiers in GRPS.

% -----Author:-----
% by David L. Jones, April-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

grps    = grps(:); % col vector
[nr,nc] = size(x);

if (nr ~= size(grps,1))
   error('# of rows in X and GRPS must be  equal !');
end

groups = unique(grps);   % unique groups
groups = groups(:);      % col vector
noGrps = size(groups,1); % # of groups 

gMean = zeros(noGrps,nc); % preallocate

for i=1:noGrps
   idx = find(grps == groups(i));      % get indices for members group i
   gRel(i,:) = sum(x(idx,:))./sum(x);  % proportion of this group to total
end

