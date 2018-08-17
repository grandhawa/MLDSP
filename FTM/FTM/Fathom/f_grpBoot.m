function B = f_grpBoot(X,grp,c)
% - within-group bootstrap sampling (fixed size)
%
% USAGE: B = f_grpBoot(X,grp,c)
%
%  X   = input data (rows = observations, cols = variables)
% grp  = column vector of integers specifying group membership
% c    = resample each column separately (default = 0);
%
% B    = boostrapped version of X
%
% SEE ALSO: f_grpResample, f_boot, f_shuffle

% -----Notes:-----
% This function performs bootstrapped resampling with replacement. The
% resampling is performed separately for each group and can also optionally be
% performed separately for each column (variable) within each group. This last
% option reduces the level of duplicate observations created in the bootstrapped
% sample.

% -----References:-----
% Fu, W. J., R. J. Carroll, and S. Wang. 2005. Estimating misclassification
%   error with small samples via bootstrap cross-validation. Bioinformatics
%   21(9): 1979-1986.

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), c = 0; end % default bootstrap columns separately

grp = grp(:);    % force col vector
n   = size(X,1); % # of rows (observations)

if (n ~= numel(grp))
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
   if (c>0) % Resample each column separately:
      B(idx.G{i},:) = f_boot(X(idx.G{i},:),1);
      
   else % Resample observations:
      % Initialize each group with 3 distinct obs (Fu et al., 2005):
      boot               = f_shuffle(X(idx.G{i},:),4);  % shuffle obs of this group
      B(idx.G{i}(1:3),:) = boot(1:3,:);                 % initialize with 3 of these
      
      % Fill in remaining obs:
      if idx.n{i}>3
         B(idx.G{i}(4:end),:) = f_boot(boot,0,idx.n{i}-3); % boostrap obs
      end
   end
end

