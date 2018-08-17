function  Y = f_cellstr2cat(X)
% - convert a cell array of strings to a single categorical variable
%
% USAGE: Y = f_cellstr2cat(X)
%
% X = cell array of strings
% Y = column vector of integers specifying factor levels or group membership
% 
% SEE ALSO: f_dummy2cat, f_dummy

% -----Author:-----
% by David L. Jones, Jun-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Notes:-----
% The purpose of this code is to create a categorical variable from a cell
% array of strings representing Gregorian dates (or some other text
% labeling system). This can then be used as a grouping variable in MANOVA,
% f_dummy, f_grpMean, etc.

% -----Check input:-----
if ~iscell(X), error('X must be a cell array of strings!'); end
% ----------------------

X  = X(:);        % force column vector
nr = size(X,1);   % get # rows

uG = f_unique(X); % unique groups, unsorted
nG = numel(uG);   % # groups
Y = nan(nr,1);    % preallocate

for i=1:nG
   idx    = ismember(X,uG(i)); % index to rows of X belonging to this group
   Y(idx) = i;
end
