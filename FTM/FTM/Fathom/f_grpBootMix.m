function B = f_grpBootMix(X,grp,n)
% - within-group bootstrap sampling (variable size), with mixing among columns
%
% USAGE: f_grpBoot(X,grp,n)
%
%  X   = input data (rows = observations, cols = variables)
% grp  = column vector of integers specifying group membership
% n    = number of new observations per groups
% 
% xB    = bootstrapped version of X
% xGrp  = bootstrapped version of grp
% 
% SEE ALSO: f_boot

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
grp = grp(:);    % force col vector
nr  = size(X,1); % # of rows (observations)

if (nr ~= numel(grp))
   error('# of rows in X and GRP must be equal !');
end
% -------------------------------------

uGrp   = f_unique(grp);  % unique groups, unsorted
noGrp  = size(uGrp,1);   % # of groups

% Preallocate:
idx.G{noGrp} = NaN;
idx.n{noGrp} = NaN;
B            = repmat(NaN,size(X));

% Get indices to members of each group:
for i=1:noGrp
   idx.G{i} = find(grp == uGrp(i));
   idx.n{i} = numel(idx.G{i});
   if idx.n{i} < 3
      error('Each group must have at least 3 members!')
   end
end

% Create boostrap sample:
for i = 1:noGrp
   % Initialize each group with 3 distinct obs:
   ini = f_shuffle(1:idx.n{i});             % pick 3 random indices to this group
   ini = ini(1:3);                          % keep top 3
   B(idx.G{i}(1:3),:) = X(idx.G{i}(ini),:); % initialize first 3 obs in this group
   
   % Fill in remaining obs:
   if idx.n{i}>3
      B(idx.G{i}(4:end),:) = X(f_boot(idx.G{i},idx.n{i}-3),:); % boostrap obs
   end
end
