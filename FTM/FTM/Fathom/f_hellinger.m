function H = f_hellinger(x,dis)
% - Hellinger transform data or create a symmetric dissimilarity matrix
% 
% USAGE: H = f_hellinger(x,dis)
% 
% x   = input matrix                (rows = observations, cols = variables)
% dis = return dissimilarity matrix (default = 0)
% 
% H   = Hellinger-transformed data or dissimilarity matrix
% 
% SEE ALSO: f_dis

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: replaced f_euclid with f_dis
% Apr-2013: updated documentation

% -----References:-----
% Legendre, P. and E. D. Gallagher. 2001. Ecologically meaningful
%   transformations for ordination of species data. Oecologia 129: 271-280.

% -----Set defaults:-----
if (nargin < 2), dis = 0; end; % default no dissimilarity matrix

[nr,nc]  = size(x);
H        = zeros(nr,nc);     % preallocate
rowSum   = sum(x,2);
idx      = find(rowSum > 0); % skip these to prevent divide-by-zero error
H(idx,:) = sqrt(x(idx,:) ./ repmat(rowSum(idx),1,nc)); % Hellinger-transformed

% Optionally return a Hellinger dissimilarity matrix:
if (dis>0)
   H = f_dis(H,'euc');
end

