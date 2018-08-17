function [gMean,gStdv,gSE,grp] = f_grpMean(x,grps)
% - returns the mean/stdv/se of X (column-wise) separately for groups
%
% USAGE: [gMean,gStdv,gSE,grp] = f_grpMean(x,grps);
%
% x    = column vector or matrix of input data
% grps = column vector of integers specifying group memberhip
%
% gMean = mean of each group (row-wise) for each variable in X (col-wise)
% gStdv = stdv of each group (row-wise) for each variable in X (col-wise)
% gSE   =   SE of each group (row-wise) for each variable in X (col-wise)
% grp   = corresponding group designation
%
% SEE ALSO: f_grpRel, f_grpOutlier, f_grpSize, f_plotError

% -----Notes:-----
% This function is used to calculate the mean of each variable in matrix X
% (column-wise) separately for each group specified by GRPS. The number of
% rows in gMean will equal the number of unique groups specified in GRPS.

% -----Author:-----
% by David L. Jones, April-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2006: added support for gStdv (i.e., to support f_plotError, etc.)
% Jun-2007: added support for gSE
% Jan-2008: updated documentation
% Oct-2010: now maintains original sort order of GRPS; now outputs GRP 

grps    = grps(:); % col vector
[nr,nc] = size(x);

if (nr ~= size(grps,1))
   error('# of rows in X and GRPS must be  equal !');
end

groups = f_unique(grps); % unique groups, unsorted
groups = groups(:);      % col vector
noGrps = size(groups,1); % # of groups 

gMean = zeros(noGrps,nc); % preallocate
gStdv = zeros(noGrps,nc);
gSE   = zeros(noGrps,nc);

for i=1:noGrps
   idx = find(grps == groups(i));   % get indices for members group i
   gMean(i,:) = mean(x(idx,:));     % MEAN for members of this group
	gStdv(i,:) = std(x(idx,:));      % STDV for members of this group
   gSE(i,:)   = f_stdErr(x(idx,:)); % stdErr for members of this group
end

% Rename for output:
grp = groups;