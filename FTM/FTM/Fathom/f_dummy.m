function  X = f_dummy(Y,trim)
% - dummy coding of categorical variables
%
% USAGE: X = f_dummy(Y,trim)
%
% Y    = column vector of integers specifying factor levels or group membership;
% trim = trim last column to avoid a singular matrix (default = 1)
% 
% SEE ALSO: f_dummy2cat, f_xMatrix, f_helmert, f_designMatrix, f_modelMatrix


% -----Author:-----
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), trim = 1; end % default to trimmed codes

if size(Y,2) > 1
   error('Y must be a column vector!');
end
% -------------------------------------

Y  = f_recode(Y); % force as consecutive integers
nr = size(Y,1);   % # rows
G  = unique(Y);   % groups
ng = size(G,1);   % # unique groups

X = zeros(nr,ng); % preallocate
for i = 1:ng
   X(Y == G(i),i) = 1; % obs belongs to this group
end

if trim>0 % trim last column to avoid a singular matrix
   X = X(:,1:end-1);
end
