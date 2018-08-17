function S = f_dis2sim(D,f,e)
% - convert square symmetric distance matrix to a similarity matrix
%
% USAGE: S = f_dis2sim(D,f,e);
%
% D = symmetric distance matrix
% f = weighting function: linear (= 1), concave-down (= 2), concave-up (= 3)
% e = exponent to use for f = 2 or 3
%
% S = similarity matrix with values ranging from 0-1
%
% SEE ALSO: f_eigenMaps

% -----Notes:-----
% Similarity matrices can be interpreted in terms of (1) a weighted graph, where
% nonzero values indicate the intensity of the connection among sites and (2) as
% a spatial weighting matrix specifying the potential strength of the
% interaction among sites. For construction of eigenvector maps, the weighting
% matrix 'S' is multiplied by a binary connectivity matrix 'B' to produce a
% spatial weighting matrix 'W'.
%
% Note that f=2, e=1, the same weighting function is produced as f=1. This which
% creates similarities that vary linearly with geographic distance.
% 
% 
% To create PCNM's use:
% A    = f_dis2sim(mst.tDis,2,2); % eq. 3 from Dray et al., 2006
% PCNM = f_eigenMaps(mst,A,0,0);

% -----References:-----
% Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial modelling: a
%   comprehensive framework for principal coordinate analysis of neighbour
%   matrices (PCNM). Ecological Modelling 196: 483-493.

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), f = 1; end; % default use linear weighting function
if ((nargin < 3) && (f>1))
   error('You must specify a value of E when F is 2 or 3!');
end

if (f_issymdis(D) == 0)
   error('Requires square symmetric distance matrix');
end

% Convert to similarities:
switch f
  
   case 1 % linear:
      S = 1 - (D/max(D(:)));

   case 2 % concave-down:
      S = 1 - (D/max(D(:))).^e;

   case 3 % concave-up:
      S = 1/(D.^e);
  
   otherwise
      error('F must be 1, 2, 3!')
end
