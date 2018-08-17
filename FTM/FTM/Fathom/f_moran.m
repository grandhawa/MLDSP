function [MC,p] = f_moran(Y,W,iter)
% - Moran's coefficient of spatial autocorrelation
% 
% USAGE: [MC,p] = f_moran(Y,W,iter)
%
% Y    = response variable
% W    = spatial weighting matrix created by f_eigenMaps
% iter = # iterations for permutation test
% 
% MC = Moran's coefficient
% p  = permutation-based significance
% 
% SEE ALSO: f_eigenMaps, f_eigenMapsStepwise, f_variogram

% -----References:-----
% Bivand, R. et al. 2008. Spatial dependence: weighting schemes, statistics and
%   models. The spdep package for R. Version 0.4-20. Available from
%   http://www.r-project.org
% Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial-modelling: a
%   comprehensive framework for principal coordinate analysis of neighbor
%   matrices (PCNM). Ecological Modelling 196: 483-493.
% Griffith, D. A. and P. R. Peres-Neto. 2006. Spatial modeling in ecology: the
%   flexibility of eigenfunction spatial analysis. Ecology 87(10): 2603-2613.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed. Elsevier
%   Science BV, Amsterdam. 

% -----Notes:-----
% This function tests either positive or negative spatial autocorrelation. When
% creating Moran's eigenvector maps (MEM's), you need to partition them
% according to their eigenvalues and run separate tests for the eigenvectors
% assocated with positive eigenvalues (positive autocorrelation) and those
% associated with negative eigenvalues (= negative autocorrelation).
% 
% Moran's coefficient usually ranges from -1 to +1 with positive values
% corresponding to positive autocorrelation and negative values corresponding to
% negative autocorrelation. Moran's I is more of a correlation coefficient while
% Geary's c is a distance-type function.

% -----Author:-----
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 3), iter = 0; end % default no permutation test

if (size(Y,1)~=size(W,1)), error('Size mismatch between X & W!'); end
% -------------------------------------

Y   = f_center(Y);
n   = size(Y,1);
uno = ones(n,1);
gW  = sum(sum(W)); % global sum of spatial weights (moran.R)
SSt = trace(Y'*Y); % Sum of Squares total

% Loosely based on 'moran.R' in the SPDEP package:
MC = ((n/gW) * (uno'*(Y*Y'.*W)*uno))/SSt;

% After Peres-Neto's StepWiseMoran.m (Y can be a vector or matrix):
% MC = n*(ones(n,1)'*(Y*Y'.*W)*ones(n,1))/sum(sum(W))/sum(sum(Y.*Y));

% Equation q. 5 of Dray et al., 2006):
% I   = eye(n,n);
% G   = (I-(1/n)*uno*uno')*W*(I-(1/n)*uno*uno'); % Gower's centered matrix
% MC  = (n/(uno'*W*uno))*((Y'*G*Y)/(Y'*(I-(1/n)*uno*uno')*Y));

% Griffith & Peres-Neto, 2006: p.2606):
% MC = (n/(uno'*W*uno)) * Y' * G * Y; 

% Appendix A of Griffith & Peres-Neto, 2006:
% Y  = f_center(Y);
% MC = (Y'*W*Y)/(Y'*Y);


if (iter>0)
   
   MCperm = zeros(iter-1,1);    % preallocate result array

   for i = 1:(iter-1)           % observed value is considered a permutation
      Yperm     = f_shuffle(Y); % permute rows of residuals
      MCperm(i) = f_moran(Yperm,W,0);
   end
   
   if (MC<=0) % handle negative MC's separately
      j  = find(MCperm <= MC);
   else
      j  = find(MCperm >= MC); % get permuted stats >= to observed statistic
   end

   p  = (length(j)+1)./(iter); % count values & convert to probability

else
   p  = NaN;
end

