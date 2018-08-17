function [centroid,SE,mDis] = f_centroid(X,grps,method)
% - returns coordinates of the centroid of X, optionally partitioned into groups
%
% USAGE: [centroid,SE,mDis] = f_centroid(X,grps,method);
%
% X      = n-dimensional coordinates (rows = observations,
%          (cols = dimensions
% grps   = optional vector of integers specifying group membership
%          e.g., grps = [1 1 2 2 3 3 3];
% method = use mean (= 1, default) or median (= 2);
%
% centroid = spatial mean or median
% SE       = corresponding standard error
% mDis     = maximum distance-to-centroid (group-wise)

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 08-Dec-2002: made grouping vector optional, error checking
% Nov-2010: replaced 'unique' with 'f_unique'
% Dec-2012: added support for SE
% Jul-2013: added support for mDis; code generalized to handle one or more
%           groups

% -----Check input & set defaults:-----
if (nargin < 2), grps   = []; end % default empty grouping vector
if (nargin < 3), method = 1;  end % default use mean

% Create default vector specifying a single group:
if isempty(grps)
   grps = ones(size(x,1),1);
end

if (size(X,1)<2)
   error('You need at least 2 points to compute the centroid!');
end
% -------------------------------------

grps   = grps(:);        % make sure it's a column vector
uGrps  = f_unique(grps); % f_unique grps, unsorted
noGrps = length(uGrps);  % number of unique groups

centroid(noGrps,size(X,2)) = 0; % preallocate results array
SE(noGrps,size(X,2))       = 0;
mDis                       = nan(noGrps,1);

for i=1:noGrps
   idx = find(grps==uGrps(i)); % get indices of rows to extract
   
   switch method
      case 1
         centroid(i,:) = mean(X(idx,:),1);         % compute centroid for this group
         SE(i,:)       = f_stdErr(X(idx,:));       % standard error of the mean
      case 2
         centroid(i,:) = median(X(idx,:),1);       % compute centroid for this group
         SE(i,:)       = f_stdErr(X(idx,:))*1.253; % standard error of the median
      otherwise
         error('Invalid METHOD!');
   end
   
   % Get maximum distance-to-centroid within each group:
   mDis(i) = max( sum( (repmat(centroid(i,:),numel(idx),1) - X(idx,:)).^2 ,2) ); % sum-of-squared distances
end
