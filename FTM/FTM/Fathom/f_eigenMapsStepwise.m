function [glob,cond,model] = f_eigenMapsStepwise(Y,X,mem,iter,alpha,verb)
% - stepwise selection of Moran's Eigenvector Maps (MEM's)
%
% USAGE: [glob,cond,model] = f_eigenMapsStepwise(Y,X,mem,iter,alpha,verb);
%
% Y     = matrix of response variables                 (rows = obs, cols = vars)
% X     = matrix of explanatory variables
% mem   = structure of output from f_eigenMaps
% iter  = # iterations for permutation tests                       (default = 0)
% alpha = significance level                                    (default = 0.05)
% verb  = optionally send result to display                        (default = 1)
%
% glob = structure of global effects (ALL VARIABLES) with the following fields:
%  .MC = Moran's coefficient of spatial autocorrelation of the residuals
%  .p  = p-value
% 
% cond = structure of conditional effects with the following fields:
%  .MC  = Moran's coefficient
%  .p   = p-value
%  .var = variable labels
%
% model = structure of marginal effects of selected spatial model:
%  .MC  = Moran's coefficient
%  .p   = p-value
%  .var = variable labels
% 
% SEE ALSO: f_eigenMaps, f_rdaStepwise

% -----Notes:-----
% This program determines whether there is significant spatial autocorrelation
% in the residuals from the regression of Y on X. These results are returned in
% GLOB. If it is not significant, the procedure stops. If it is significant,
% this function performs a stepwise procedure to incrementally remove the
% spatial autocorrelation in the residuals. This is accomplished by sequentially
% adding one of the spatial descriptors in MEM to X, regressing Y on the new X,
% and assessing whether there still remains significant spatial autocorrelation
% in the new residuals. The spatial variables in MEM are ranked and selected
% based upon their ability to account for spatial structure in the regression 
% model; i.e., the 'best' variable in each step of the selection process is the
% one that yields the lowest Moran's coefficient in the resulting residuals.
% Addition of spatial variables to the regression model is halted when there is
% no longer any significant spatial autocorrelation in the model's residuals.

% -----References:-----
% Griffith, D. A. and P. R. Peres-Neto. 2006. Spatial modeling in ecology: the
%   flexibility of eigenfunction spatial analysis. Ecology 87(10): 2603-2613.

% -----Author:-----
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), iter  = 0;    end % no permutation test by default
if (nargin < 4), alpha = 0.05; end % default significance level
if (nargin < 5), verb  = 1;    end % send output to display by default

% Extract components of MEM:
S = mem.evects; % spatial descriptors
C = mem.W;      % connectivity matrix

% Default S labels:
xLabels = cellstr(num2str([1:size(S,2)]'))';

% Check input sizes:
n   = size(Y,1); % # of observations
if ( n ~= size(X,1) || n ~= size(S,1) || n ~= size(C,1) )
   error('X, Y, & mem.W need same # of rows')
end

% Check if connectivity matrix is symmetric:
if (f_issymdis(C,0) == 0)
   error('C must be square symmetric matrix');
end

tol_p = 0.001;   % set tolerances
ncS   = size(S,2); % # spatial descriptors

% =========================================================================
% GLOBAL TEST:
fprintf('\nTesting GLOBAL Moran''s coefficient with %d permutations...\n',iter);

% Get residuals:
result = f_rda(Y,X,[0],0,0,0,0);
Yres   = result.res; clear result;

% Check for signifiant spatial autocorrelation in residuals:
% [glob.MC,glob.p] = s_globalMoran(Yres,C,iter);

[glob.MC,glob.p] = f_moran(Yres,C,iter);


% Assess stopping criteria:
aFlag = ([glob.p - alpha] < tol_p); % alpha flag


if (aFlag) % Significant global test:
   fprintf(['\nGLOBAL TEST (p = ' num2str(glob.p)...
      ') is significant at alpha = ' num2str(alpha)]);
   fprintf('\n\nNow performing stepwise variable selection... \n');

else % Non-significant global test:
   fprintf(['\nGLOBAL TEST (p = ' num2str(glob.p)...
      ') is NOT significant at alpha = ' num2str(alpha)]);
   fprintf('\nNo further steps are required. \n');
   
   % No variable selection:
   cond  = NaN;
   model = NaN;
   return
end

% =========================================================================
% CONDITIONAL EFFECTS:
fprintf('\nPerforming CONDITIONAL TESTS with %d permutations...\n',iter);
for i = 1:ncS
   fprintf('  -> examining: %s \n',xLabels{i});
   
   % Get residuals:
   result = f_rda(Y,[X S(:,i)],[0],0,0,0,0);
   Yres   = result.res; clear result;
   
   % Check for spatial autocorrelation in residuals:
   [MC,p]      = f_moran(Yres,C,0);
   cond.MC(i)  = MC;         % Moran's Coefficient
   cond.p(i)   = p;          % p-value
   cond.var(i) = xLabels(i); % variable label
end

% Determine the single best variable:
idx = find(cond.MC == min(cond.MC));  % smallest MC wins
if (length(idx)>1), idx = idx(1); end % in case of ties, select one

% Re-do conditional test on best variable to get p-value:
result             = f_rda(Y,[X S(:,idx)],[0],0,0,0,0);
Yres               = result.res; clear result;
[null,cond.p(idx)] = f_moran(Yres,C,iter); % replace with permuted p-value

% =========================================================================
% MARGINAL EFFECTS:
fprintf('\nPerforming MARGINAL TESTS with %d permutations...\n',iter);
noVars = ncS; % keep a copy for fprintf's
fprintf('  -> selecting variable 1 of %d \n',noVars);

% Initialize model with best variable from conditional test:
model.MC(1)  = cond.MC(idx);
model.p(1)   = cond.p(idx);
model.var(1) = cond.var(idx);

X            = [X S(:,idx)];  % add best variable to X
S(:,idx)     = [];            % remove best variable from pool
xLabels(idx) = [];            % label removed from pool

for j = 2:ncS                 % 1st variable selected from Conditional tests
   fprintf('  -> selecting variable %d of %d \n',j,noVars);
   % Marginal Tests:
   nc_s = size(S,2); % update # of remaining spatial variables
   for i = 1:nc_s

      % Get residuals:
      result = f_rda(Y,[X S(:,i)],[0],0,0,0,0);
      Yres   = result.res; clear result;

      % Check for spatial autocorrelation in residuals:
      [MC,p]      = f_moran(Yres,C,0);
      marg.MC(i)  = MC;         % Moran's Coefficient
      marg.p(i)   = p;          % p-value
      marg.var(i) = xLabels(i); % variable label
      
   end
   
   % Determine the single best variable:
   idx = find(marg.MC == min(marg.MC));   % smallest MC wins
   if (length(idx)>1), idx = idx(1); end; % in case of ties, select one

   % Re-do marginal test on best variable to get p-value:
   result             = f_rda(Y,[X S(:,idx)],[0],0,0,0,0);
   Yres               = result.res; clear result;
   [null,marg.p(idx)] = f_moran(Yres,C,iter); % replace with permuted p-value
   
   % Assess stopping criteria:
   aFlag = ([marg.p(idx) - alpha] < tol_p); % alpha flag
   
   if (aFlag)

      % Add this variable to model:
      model.MC(j)  = marg.MC(idx); % Moran's Coefficient
      model.p(j)   = marg.p(idx);  % p-value
      model.var(j) = xLabels(idx); % variable label
      
      % Prepare for next iteration:
      X            = [X S(:,idx)]; % add best variable to X
      S(:,idx)     = [];           % remove best variable from pool
      xLabels(idx) = [];           % label removed from pool
      clear marg;                  % clean up for next iteration
 
   else % Don't add any more variables to model
      break
   end
end
% =========================================================================


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('Stepwise selection of MORAN''S EIGENVECTORS:\n');
   fprintf('--------------------------------------------------\n');

   fprintf('Global Test: (all variables included)\n');
   headerCell = {'MC','p'};    % set up column labels
   resultCell = f_struct2flat(glob);        % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('\nUsed %d permutations of residuals \n',iter);
   fprintf('--------------------------------------------------\n\n');

   fprintf('Conditional Tests: (each variable separately)\n');
   headerCell = {'MC','p', 'Variable'};     % set up column labels
   resultCell = f_struct2flat(cond);        % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('\nUsed %d permutations of residuals \n',iter);
   fprintf('--------------------------------------------------\n\n');

   fprintf('Marginal Tests: (sequential variable addition)\n');
   headerCell = {'MC','p','Variable'};      % set up column labels
   resultCell = f_struct2flat(model);       % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('\nUsed %d permutations of residuals \n',iter);

   fprintf('(Only the selected spatial variables are shown)\n');
   fprintf('--------------------------------------------------\n\n');

   % Report if stopping criterion was met:
   if (~aFlag)
      fprintf('Variable addition halted due to: ALPHA LEVEL. \n');
   else
      fprintf('Stopping criterion was not met. \n');
   end

end
% ---------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                SUBFUNCTIONS:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MC,p] = s_globalMoran(Y,W,iter)
% - calculate global Moran's Coefficient
%
% Y    = response variable (e.g., residuals from regressing Y on X)
% W    = spatial weighting matrix
% iter = # iterations for permutation test
% 
% result = structure with the following fields:
%   .MC = global Moran's Coefficient
%   .p  = p-value

% Global Moran's coefficient:
MC = (Y' * W * Y) / (Y' * Y); % Griffith & Peres-Neto, 2006: Appendix A):


% Permutation test:
if (iter>0)

   MCperm = zeros(iter-1,1);    % preallocate result array

   for i = 1:(iter-1)           % observed value is considered a permutation
      Yperm     = f_shuffle(Y); % permute rows of residuals
      MCperm(i) = s_globalMoran(Yperm,W,0);
   end

   if MC<0 % handle negative MC's separately
      j  = find(MCperm <= MC);
   else
      j  = find(MCperm >= MC); % get permuted stats >= to observed statistic
   end

   p  = (length(j)+1)./(iter); % count values & convert to probability

else
   p  = NaN;
end


