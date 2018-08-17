function [n,uGrp] = f_grpSize(grp)
% - get the number of observations in each category for a grouping variable
%
% USAGE: function [n,uG] = f_grpSize(grp)
%
% grp = column vector of whole numbers specifying group memberhip
% 
% n    = # observations in each group
% uGrp = corresponding group
%
% SEE ALSO: f_grpMean

% -----Notes:-----
% This function determines the number of observations represented in each
% category of a grouping variable used for MANOVA, ANCOVA, Discriminant
% Analysis, etc. I allows you quick way to determine whether your
% experimental design is balanced or not.

% -----Author:-----
% by David L. Jones, July-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

grp  = grp(:);        % force col vector
uGrp = f_unique(grp); % unique groups, unsorted
nGrp = numel(uGrp);   % # of groups 
n    = zeros(nGrp,1); % preallocate

for i=1:nGrp
   n(i) = sum(grp==uGrp(i)); % # observations in this group
end
