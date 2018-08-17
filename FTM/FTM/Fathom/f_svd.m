function [U,evals,V] = f_svd(x)
% - singular value decomposition of a matrix
%
% USAGE: [U,evals,V] = f_svd(x);
%
% x     = input matrix
% U     = left eigenvectors
% evals = sorted eigenvalues
% V     = right eigenvectors
%
% SEE ALSO: f_eig

% -----Notes:-----
% SVD does not return negative eigenvalues as EIG does.

% -----Author(s):-----
% by David L. Jones, Sep-2003
% 
% This file is part of the FATHOM Toolbox for Matlab and 
% is released under the GNU General Public License, version 2.

% Mar-2011: slight overhaul, now outputs V, and trims U,V

% Singular Value Decomposition:
[U,D,V] = svd(x); 
evals   = diag(D);

% Trim axes:
s = numel(evals);
U = U(:,1:s);
V = V(:,1:s);


