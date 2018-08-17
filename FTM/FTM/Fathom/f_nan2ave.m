function Y = f_nan2ave(X,grp)
% - replace missing values (NaN's) with average value (columnwise) by group
%
% USAGE: Y = f_nan2ave(X,grp)
% 
% X   = input matrix (rows are observations, columns are variables)
% grp = optional grouping vector
% 
% Y   = output with NaN's replaced by group means

% -----Author:-----
% by David L. Jones, Feb-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Sep-2014: now supports multiple columns and an optional grouping vector

% -----Check input & set defaults:-----
if (nargin <  2), grp = ones(size(X,1),1); end % default 1 group

% Check grouping vector:
grp = grp(:);
if size(X,1) ~= numel(grp)
   error('X & GRP must have the same # of rows!')
end
% -------------------------------------

Y    = X;             % make a copy
uGrp = f_unique(grp); % unique groups, unsorted
nGrp = numel(uGrp);   % # of groups
nc   = size(X,2);     % # columns  

for i=1:nGrp % repeat for each group
   
   % Get index to rows of group i:
   idxG = find(grp==uGrp(i)); 
      
   for j = 1:nc % repeat for each column
      % Get index to rows of group i of column j that are NaN's:
      idxN = find( (grp==uGrp(i)) & isnan(X(:,j)));
      
      % Replace NaN's with group mean:
      Y(idxN,j) = nanmean(X(idxG,j));
   end
end
