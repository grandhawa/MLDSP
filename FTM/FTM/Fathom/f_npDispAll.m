function result = f_npDispAll(yDis,grps,sm,verb)
% - similar to f_npDisp, but returns distance to all centroids
%
% USAGE: result = f_npDispAll(yDis,grps,sm,verb)
%
% yDis = square symmetric dissimilarity matrix derived from response variables
% grps = column vector of integers specifying group membership for objects in yDis
% sm   = use spatial median instead of centroid              (default = 1)
% verb = optionally send results to display                  (default = 1)
%
% result = structure with the following fields:
%  .z       = residuals (= distance each observation is from its group's spatial
%             median or centroid)
%  .grps    = grouping vector
%  .gLabels = group labels
%  .scores  = scaled PCoA eigenvectors
%  .evals   = eigenvalues
%  .C       = coordinates of each group's (row-wise) spatial median or centroid
%             in REAL space
%  .gMean   = mean of the residuals within each group
%  .type    = 'spatial median' or 'centroid'
% 
% SEE ALSO: f_npDisp

% -----Notes:-----
% Currently this code only supports 1 grouping variable
% 
% The residuals obtained from this function can subsequently be used to create
% a Eulcidean distance matrix that, along with the grouping vector, can be used
% as input to 'f_npManova' to perform a permutation-based significance test of
% differences in multivariate dispersion among the groups.

% -----References:-----
% Anderson, M. J. 2006. Distance-based tests for homogeneity of multivariate
%   dispersions. Biometrics 62: 245-253
% Anderson, M. J., K. E. Ellingsen, and B. H. McArdle. 2006. Multivariate
%   dispersion as a measure of beta diversity. Ecology Letters 9(6): 683-693

% -----Author:-----
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2010: unique replaced with f_unique
% Jan-2011: uGrps, noGrps calculated first

% -----Check input & set defaults:-----
if (nargin < 3), sm   = 1; end; % default use spatial median
if (nargin < 4), verb = 1; end; % send output to display by default

n = size(yDis,1); % # observations

if n ~= size(grps,1), error('yDis & GRPS need same # of rows'); end;

if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix');
end

if (~size(grps,2)==1)
   error('Only 1 grouping vector is currently supported!')
end
% -------------------------------------

% Set up groups:
uGrps  = f_unique(grps); % unique grps, unsorted
noGrps = length(uGrps);  % number of unique groups

% Get Principal Coordinates (scaled), keep negative eigenvalues:
[U,evals] = f_pcoa(yDis,0,1,1);

% Get # of negative eigenvalues:
neg = sum(find(evals<0));

z = repmat(NaN,n,1); % initialize
if (neg>0) %--Negative Eigenvalues:--
   
   % Extract PC's associated with negative eigenvalues:
   idx          = find((evals<0)==1); % get index to negative eigenvalues
   U_neg        = U(:,idx);           % extract PC's with neg eigenvalues
   U_pos        = U;                  % make a copy
   U_pos(:,idx) = [];                 % keep PC's with pos eigenvalues
   
   C_pos = repmat(NaN,noGrps,size(U_pos,2)); % preallocate results array
   C_neg = repmat(NaN,noGrps,size(U_neg,2)); % preallocate results array
   
   for i = 1:noGrps
      idx = find(grps==uGrps(i)); % get indices of rows to extract
      
      % Get spatial median/centroids of each group:
      if (sm>0)
         C_pos(i,:) = median(U_pos(idx,:)); % REAL space
         C_neg(i,:) = median(U_neg(idx,:)); % IMAGINARY space
      else
         C_pos(i,:) = mean(U_pos(idx,:)); % REAL space
         C_neg(i,:) = mean(U_neg(idx,:)); % IMAGINARY space
      end
      
      % Get distance of each observation to centroid (residuals):
      r_pos = f_dis([C_pos(i,:);U_pos(idx,:)],'euc'); % REAL SPACE
      r_pos = r_pos(2:end,1);
      
      r_neg = f_dis([C_neg(i,:);U_neg(idx,:)],'euc'); % IMAGINARY SPACE
      r_neg = r_neg(2:end,1);
      
      % Subtract IMAGINARY distance from REAL distance:
      z(idx) = sqrt(abs(r_pos.^2 - r_neg.^2)); % ABS() avoids imaginary distances when r_neg > r
   end
   
else %--No Negative Eigenvalues:--
   % Make a copy:
   U_pos = U;
   
   C_pos = repmat(NaN,noGrps,size(U_pos,2)); % preallocate results array
      
   for i = 1:noGrps
      idx = find(grps==uGrps(i)); % get indices of rows to extract
      
      % Get spatial median/centroids of each group:
      if (sm>0)
         C_pos(i,:) = median(U_pos(idx,:)); % REAL space
      else
         C_pos(i,:) = mean(U_pos(idx,:));   % REAL space
      end
      
      % Get distance of each observation to centroid (residuals):
      r_pos = f_dis([C_pos(i,:);U_pos(idx,:)],'euc'); % REAL SPACE
      r_pos = r_pos(2:end,1);
      
      z(idx) = r_pos;
   end
end

% Group means:
gMean = f_grpMean(z,grps);

% Wrap results up into a structure:
result.z       = z;
result.grps    = grps;
result.gLabels = cellstr(num2str(f_unique(grps))); % group labels, unsorted
result.scores  = U;
result.evals   = evals;
result.C       = C_pos;
result.gMean   = gMean;
if (sm>0)
   result.type = 'spatial median';
else
   result.type = 'centroid';
end


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('      Homogeneity of Multivariate Dispersions:');
   fprintf('\n==================================================\n\n');
   fprintf('# Pos Eigenvalues = %d\n', sum(evals>=0));
   fprintf('# Neg Eigenvalues = %d\n\n', sum(evals <0));
   txt = ['Average distance to ' result.type ':\n'];
   fprintf(txt);
   for i = 1:noGrps
      fprintf('Group %d = %2.4f\n',uGrps(i),gMean(i));
   end
end
% ---------------------------------
