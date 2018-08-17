function result = f_covPool(X,grp)
% - covariance matrix pooled across groups
%
% USAGE: result = f_covPool(X,grp)
%
% X   = matrix of input data (row = observations, cols = variables)
% grp = column vector of integers specifying group membership
%
% result = pooled covariance matrix
%
% SEE ALSO: cov

% -----References:-----
% http://www.mathkb.com/Uwe/Forum.aspx/matlab/16779/pooled-variance-covariance-matrix

% -----Author:-----
% by David L. Jones, Jan-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (size(X,1) ~= size(grp,1)), error('Size mismatch b/n X & xGrp'); end
% -------------------------------------

% Pooled covariance:
result = cov(f_center(X,grp));
   
