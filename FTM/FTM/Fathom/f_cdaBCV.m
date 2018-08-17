function f_cdaBCV(X,grp,method,iter)
% - bootstrap cross-validation for Canonical Discriminant Analysis
%
% USAGE: f_cdaBCV(X,grp,method,iter)
%
% x   = input data (rows = observations, cols = variables)
% grp = column vector of integers specifying group membership
%
% method = method of classification:
%          1 : linear
%          2 : quadratic
%          3 : mahalanobis
%          4 : centroid     (default)
%
% iter  = # iterations for bootstrap resampling (default = 200)
%
% SEE ALSO: f_cdaCV, f_cda632, classify, f_grpBoot

% -----Notes:-----
% This function is used to diagnose a Canonical Discriminant Analysis using
% a new bootstrap method of cross-validation, specifically for small sample
% sizes.

% -----References:-----
% Fu, W. J., R. J. Carroll, and S. Wang. 2005. Estimating misclassification
%   error with small samples via bootstrap cross-validation. Bioinformatics
%   21(9): 1979-1986.

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), method = 4;   end % 'centroid' method by default
if (nargin < 4), iter   = 200; end % default 200 bootstraps

switch method
   case 1, type = 'linear';
   case 2, type = 'quadratic';
   case 3, type = 'mahalanobis';
   case 4, type = 'centroid';
   otherwise, error('Classification method must be 1, 2, 3, or 4!');
end

grp = grp(:);    % force col vector
n   = size(X,1); % # of rows (observations)

if (n ~= numel(grp))
   error('# of rows in X and GRP must be equal !');
end

if (iter<50), error('ITER should be from 50-200!'); end
% -------------------------------------

uGrp   = f_unique(grp);  % unique groups, unsorted
noGrp  = size(uGrp,1);   % # of groups

% Preallocate:
err.tot      = repmat(NaN,iter,1);
err.grp      = repmat(NaN,noGrp,iter);

% Repeat for each bootstrap sample:
for j = 1:iter
   
   % Create boostrap sample:
   B = f_grpBoot(X,grp,0);
   
   % Leave-One-Out Cross-Validation of Bootstrapped data:
   errB = f_cdaCV(B,grp,method,0);
      
   % Collect classification error rates for each iteration:
   err.tot(j)         = errB.tot;
   err.grp(1:noGrp,j) = errB.grp;
end

% Average results:
err.tot = mean(err.tot);
err.grp = mean(err.grp,2);

% -----Send output to display:-----

fprintf('\n==================================================\n');
fprintf('        F_CDA BOOTSTRAP CROSS-VALIDATION\n'           );
fprintf('            Classification Success: \n'               );
fprintf('--------------------------------------------------\n' );

fprintf('Group        Corrrect  \n');
for j=1:noGrp
   fprintf('%s %d %s %10.1f %s \n',['  '],uGrp(j),['     '],[1-err.grp(j)]*100,['%']);
end

fprintf('\n\n');
fprintf('Total Correct  = %4.2f %s \n',(1-err.tot)*100,['%']);
fprintf('Total Error    = % 4.2f %s \n',err.tot*100,['%']);
fprintf('Class. method  = %s \n',type);
if method<4, fprintf('Prior prob     = group size \n'); end
fprintf('\n==================================================\n\n');



