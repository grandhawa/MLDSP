function idx = f_balance(grps)
% - get indices of elements to create a balanced design
%
% USAGE idx = f_balance(grps);
%
% grps = vector of integers specifying group membership
% idx  = index of elements for a balanced design
%
% SEE ALSO: f_randSub

% -----Notes:-----
% The purpose of this function is to assist the user in creating a balanced
% data matrix, in terms of group size, for subsequent use in MANOVA, CDA,
% etc. Given a matrix of response data (X with rows = observations, columns =
% variables) and a vector of integers specifying group membership (GRPS),
% this function will return a vector of indices (IDX) that specifies a
% subset of X and GRPS that all have equal group sizes (= to the size of the
% smallest group).
%
% The variable IDX can be used in subsequent analyses requiring a balanced
% designed via X(idx,:) and grps(idx).

% -----Author:-----
% by David L. Jones, Aug-2004
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% To do: allow user to specify min group size, and make sure it is <= n.

uGrps   = unique(grps);  % unique groups
noGrps  = length(uGrps); % # unique groups   
gSize   = [];            %initialize variable

for i = 1:noGrps
   [gRows{i},ignore] = find(grps==uGrps(i)); % get row indices for each group
   gSize = [gSize; size(gRows{i},1)];        % get size of each group
end

n   = min(gSize); % get size of smallest group
idx = [];         % initialize

% Randomly subset, so all groups are same size as smallest:
for i = 1:noGrps
   sub = f_randSub(gRows{i},n);
   idx = [idx; sub];
end
