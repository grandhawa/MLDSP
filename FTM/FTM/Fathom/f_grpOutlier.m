function idx = f_grpOutlier(X,grp,tol)
% - logical index identifying outliers, separately for each group
%
% USAGE: idx = f_grpOutlier(X,grp,tol);
%
% X   = input data (rows = obs, cols = variables)
% grp = column vector of integers specifying group memberhip
%        (default all belong to the same group)
% tol = tolerance for detecting an outlier (default = 10)
%
% idx = logical index (0 = outlier, 1 = NOT an outlier)
%
% SEE ALSO: f_outlier

% -----Notes:-----
% Logical indices are returned with '1' indicating an observation is NOT an
% outlier, so the 'non-outlier' data can be used simply by calling the index.
% For example:
% > out   = f_grpOutlier(X,grp)
% > gMean = f_grpMean(X(out,:),grp(out))
% 
% Breiman & Cutler (2003) defined an outlier as an observation that only
% has weak similarities to the other members of its class.

% -----References:-----
% Breiman, L., and A. Cutler. 2003. Manual on setting up, using, and
%   understanding Random Forests v4.0. Technical Report. 
%   ftp://ftp.stat.berkeley.edu/pub/users/breiman/Using_random_forests_v4.0.pdf

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2012: updated documentation

% -----Set defaults and check input:-----
if (size(X,1)==1), X = X(:);               end % col vector
if (nargin < 2), grp =  ones(size(X,1),1); end % default only 1 class
if (nargin < 3), tol = 10;                 end % default tolerance of 10

grp  = grp(:); % col vector
[nr] = size(X,1);

if (nr ~= size(grp,1))
   error('# of rows in X and GRPS must be  equal !');
end
% ----------------------

dis  = f_dis(X,'euc');       % Euclidean distance matrix
dis  = dis/max(dis(:));      % recode so values range from 0-1
sim  = 1-dis;                % convert to similarities
idx  = ones(size(grp));      % preallocate
out  = f_outlier(sim,grp,0); % calculate outlier measures

% Identify values > tolerance as outliers:
idx(out>=tol) = 0;

% Force output as logical
idx = logical(idx);
