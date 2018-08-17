function [model,idx,result] = f_selectW(Y,N,A,neg)
% - AIC-based selection of optimal spatial weighting matrices
% 
% USAGE: [model,idx] = f_selectW(Y,N,A,neg);
% 
% Y   = response variable
% N   = cell array of neighbor graphs created by f_delaunay, f_dnn, f_gabriel,
%       f_mst, or f_relNeigh
% A   = cell array of corresponding weighting matrices created by f_dis2sim
% neg = keep negative eigenvalues  (default = 0)
% 
% model = structure of marginal effects of the optimal model:
%  .RSS   = residual sum-of-squares
%  .R2    = CUMULATIVE fraction of total variance explained
%  .R2adj = CUMULATIVE adjusted R2
%  .AIC   = corrected AIC (or BIC)
%  .var   = variable labels
%  .idx   = index to the selected variables
%  .null  = value of AIC for a null model (i.e., X is only an intercept term)
%  .X     = matrix of explanatory variables of the selected model
% 
% idx = index to optimal model
% 
% SEE ALSO: f_eigenMaps

% -----Author:-----
% by David L. Jones, Mar-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% May-2013: updated documentation

% -----Check input & set defaults:-----
if (nargin < 3), neg  = 0; end; % default don't keep negative eigenvalues

if (~iscell(N)) || (~iscell(A))
   error('N and A must be cell arrays!')
end

if size(N,2) ~= size(A,2)
   error('N and A must be the same size!')
end
% ----------------------

n         = size(N,2);  % get number of competing models to evaluate
AIC       = zeros(n,1); % preallocate
MEM{n}    = NaN; 
result{n} = NaN;

% Cycle through spatial models:
for i = 1:n
   % Show progress:
   fprintf('Examining model %d of %d \n',i,n);
   
   % Create Eigenvector Maps, no iterations, KEEP neg eigenvalues:
   MEM{i} = f_eigenMaps(N{i},A{i},0,neg);

   % AIC-based stepwise selection of MEM's:
   result{i} = f_rdaAIC(Y,MEM{i}.evects,0,0);
   
   % Tabulate results:
   AIC(i) = result{i}.minAIC;
end
 
% Determine the optimal model:
idx = find(AIC == min(AIC));
idx = idx(1); % in case of ties, select one
 
% Extract optimal model:
model = result{idx};
