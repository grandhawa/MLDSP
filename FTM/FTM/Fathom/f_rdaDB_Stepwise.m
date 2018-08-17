function [glob,cond,model] = f_rdaDB_Stepwise(yDis,ncY,X,iter,verb,alpha,stop,grp,xLabels)
% - stepwise selection of explanatory variables in db-RDA based on F-stat
%
% USAGE: [glob,cond,model] = f_rdaStepwise(yDis,ncY,X,iter,verb,alpha,stop,grp,xLabels)
%
% yDis  = square symmetric distance matrix derived from response variables
% ncY   = # columns of original transformed data used to derive yDis
% X     = matrix of explanatory variables
% iter  = # iterations for global/marginal permutation tests;      (default = 0)
%         input 2 values to speed up conditional tests;    e.g., iter = [1000 0]
% verb  = optionally send result to display                        (default = 1)
% alpha = significance level                                    (default = 0.05)
% stop  = enforce stopping criteria                                (default = 1)
% grp   = row vector of integers specifying which columns of X belong to the
%         same variable                                            (default = 0)
%         e.g., if cols 3-5 are dummy codes for a variable, grp = [1 2 3 3 3]
%            
% xLabels = cell array of X labels                        (default = autocreate)
%           e.g., xLabels = {'temp' 'sal' 'depth'}
%
% glob = structure off global effects (ALL VARIABLES) with the following fields:
%  .F     = F-stat
%  .p     = p-value
%  .R2    = fraction of TOTAL variance explained
%  .R2adj = adjusted R2
% 
% cond = structure of conditional effects with the following fields:
%  .F     = F-stat
%  .p     = p-value
%  .R2    = fraction of TOTAL variance explained
%  .R2adj = adjusted R2
%  .var   = variable labels
%
% model = structure of marginal effects of selected model:
%  .F     = partial F-stat
%  .p     = p-value
%  .R2    = fraction of (PARTIAL) variance explained
%           (after effects of previous terms in model are removed)
%  .R2adj = adjusted R2
%  .var   = variable labels
%
% SEE ALSO: f_rdaDB, f_rdaDB_AIC, f_rdaStepwise

% -----Notes:-----
% This function is used to perform selection of explanatory variables in
% db-RDA via a stepwise, forward addition procedure. Variables are
% sequentially added by selecting the one that yields the largest partial
% F-statistic. Following the recommendations of Blanchet et al., three
% stopping criteria for variable addition are used: 1) global test of
% signifance, 2) permutation-based p-value of the variable selected in each
% step, and 3) the corresponding R^2adj. Note that automatic model building
% techniques are used to provide a parsimonious, 'optimal' subset of
% available variables, but are not guaranteed to find the 'best' subset.
%
% CONDITIONAL tests examine the independent effect of each explanatory variable
% on the response matrix yDis. A second value input in ITER will use a different
% number of permutations for the conditional tests vs. the global/marginal
% tests. This is primarily used to speed up the variable selection procedure,
% since ITER = [1000 0] will skip the permutations in the conditional tests, but
% still use 1000 permutations for the others.
%
% MARGINAL tests examine the effect each explanatory variable has on the
% response matrix yDis after taking into account (partial out) the effects of
% explanatory variables selected previously during the stepwise variable
% selection procedure.
% 
% Use GRP to specify which variables in X that represent dummy-coded categorical
% variables that should be added to the model together as a single entity.
% 
% p-values can be adjusted for multiple tests by the Sidak Method using the
% formula: p_corrected = 1-(1-p)^(k), where p is the original probability and k
% is the number of tests. Blanchet et al. give an example of using this when a
% dataset is 'saturated' (# variables = n-1). In this case they spit the
% variables in half and performed 2 separate tests, but corrected the p-values
% using k=2.

% -----References:-----
% Blanchet, F. G., P. Legendre, & D. Borcard. 2008. Forward selection of
%   explanatory variables. Ecology 89:2623-2632.

% -----Author:-----
% by David L. Jones, July-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% modified after f_rdaStepwise to support db-RDA

% -----Set defaults & check input:-----
if (nargin < 3), iter    = 0;    end % no permutation test by default
if (nargin < 4), verb    = 1;    end % send output to display by default
if (nargin < 5), alpha   = 0.05; end % default significance level
if (nargin < 6), stop    = 1;    end % default enforce stopping criteria
if (nargin < 7), grp     = 0;    end % default treat all variables separate
if (nargin < 8), xLabels = cellstr(num2str([1:size(X,2)]'))'; end % default X labels

% Check for separate global/marginal and conditional 'iter':
iter = iter(:);
if size(iter,1)>2, error('Please specify only 1 or 2 values for ITER!'); end
if size(iter,1)<2, iter = [iter;iter]; end

% If labels are not cell arrays, try forcing them:
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end
xLabels = xLabels(:); % force row vector

% Check input:
if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix!');
end
% -------------------------------------

% -----Process groups:-----
if (sum(grp(:))==0)   % all explanatory variables are treated separate
   grp = 1:size(X,2); % autocreate groups
elseif (size(grp,1)~=1) || (size(grp,2)~=size(X,2))
   error('GRP must be a row vector with the same # columns as X!');
end

uGrp = unique(grp);      % get list of unique groups of explanatory variables
ng   = size(uGrp,2);     % get number of unique groups

g = zeros(ng,size(X,2)); % preallocate
for i = 1:ng
   g(i,:) = (grp == uGrp(i)); % find columns of X that match this group
end
g = logical(g);               % convert to logical indices
% -------------------------

% Check input sizes:
n   = size(yDis,1); % # of observations
if n ~= size(X,1), error('yDis & X need same # of rows'); end

if ng ~= size(xLabels,1)
   error('# of xLabels and explanatory variables do not match!');
end

% Set tolerances:
tol_p = 0.001;
tol_r = 0.00001;

% =========================================================================
% GLOBAL TEST:
fprintf('\nPerforming GLOBAL TEST with %d permutations...\n',iter(1));
result     = f_rdaDB(yDis,ncY,X,[0],iter(1),0,0,0);
glob.F     = result.F;
glob.p     = result.p;
glob.R2    = result.R2;
glob.R2adj = result.R2adj;
clear result;
if (glob.p > alpha) % Non-significant global test:
   fprintf(['\nGLOBAL TEST (p = ' num2str(glob.p)...
      ') is NOT significant at alpha = ' num2str(alpha)]);
   fprintf('\nNo further steps are required. \n');

   if (stop==1) % Enforce stopping rules:
      cond  = NaN;
      model = NaN;
      return;
   end
else % Significant global test:
   fprintf(['\nGLOBAL TEST (p = ' num2str(glob.p)...
      ') is significant at alpha = ' num2str(alpha)]);
   fprintf('\n\nNow performing stepwise variable selection... \n');
end


% =========================================================================
% CONDITIONAL EFFECTS:
fprintf('\nPerforming CONDITIONAL TESTS with %d permutations...\n',iter(2));
for i = 1:ng
   fprintf('  -> examining: %s \n',xLabels{i});
   result        = f_rdaDB(yDis,ncY,X(:,g(i,:)),[0],iter(2),0,0,0);
   cond.F(i)     = result.F;     % F-stat
   cond.p(i)     = result.p;     % p-value
   cond.R2(i)    = result.R2;    % fraction of TOTAL variation explained
   cond.R2adj(i) = result.R2adj; % adjusted R^2
   cond.var(i)   = xLabels(i);   % variable label
end

% Determine the single best variable:
idx = find(cond.F == max(cond.F)); % best F-ratio wins
idx = idx(1);                      % in case of ties, select one

% Repeat conditional permutation test on best variable using global/marginal iter:
best        = f_rdaDB(yDis,ncY,X(:,g(idx,:)),[0],iter(1),0,0,0);
cond.p(idx) = best.p; % replace with permuted p-value

% =========================================================================
% MARGINAL EFFECTS:
fprintf('\nPerforming MARGINAL TESTS with %d permutations...\n',iter(1));
noVars = ng; % keep a copy for fprintf's
fprintf('  -> selecting variable 1 of %d \n',noVars);

% Initialize model with best variable from conditional test:
model.F(1)      = cond.F(idx);
model.p(1)      = cond.p(idx);
model.R2(1)     = cond.R2(idx);
model.R2adj(1)  = cond.R2adj(idx);
model.seqVar(1) = cond.R2adj(idx); % for 1st variable, PARTIAL explained = TOTAL explained
model.var(1)    = xLabels(idx);

W             = X(:,g(idx,:)); % selected variable becomes covariable
X(:,g(idx,:)) = [];            % covariable removed from pool
uGrp(idx)     = [];            % variable group removed from pool
delCol        = logical(g(idx,:)==1); % find associated columns to remove
g(:,delCol)   = [];            % remove associated column
g(idx,:)      = [];            % indices to variable group removed from pool
xLabels(idx)  = [];            % label removed from pool

for j = 2:(ng)           % 1st variable selected from Conditional tests
   fprintf('  -> selecting variable %d of %d \n',j,noVars);
   % Marginal Tests:
   n_g = size(uGrp,2);      % update # of remaining variable groups
   for i = 1:n_g
      result    = f_rdaDB(yDis,ncY,X(:,g(i,:)),[W],[0],0,0,0);
      marg.F(i) = result.F;   % partial F-stat
   end

   % Select best variable from MARGINAL tests:
   idx = find(marg.F == max(marg.F));
   idx = idx(1); % in case of ties, select one

   % Re-do marginal test on best variable, but now with permutation:
   result = f_rdaDB(yDis,ncY,X(:,g(idx,:)),[W],iter(1),0,0,0);
   
   % Get cumulative fraction of TOTAL variation explained:
   result_tot = f_rdaDB(yDis,ncY,[W X(:,g(idx,:))],[0],0,0,0,0); % covariables not removed

   % Assess stopping criteria:
   aFlag = ((result.p - alpha) < tol_p);              % alpha flag
   rFlag = ((result_tot.R2adj - glob.R2adj) < tol_r); % R^2adj flag
   
   % When stop==1, both stopping flags must be true to add another variable:
   if (stop==0) || (aFlag && rFlag)

      % Add this variable to model:
      model.F(j)      = result.F;      % partial F-stat
      model.p(j)      = result.p;      % p-value
      model.R2(j)     = result.R2;     % fraction of PARTIAL variation explained
      model.R2adj(j)  = result.R2adj;  % adjusted R^2
      model.seqVar(j) = result_tot.R2adj; % fraction of TOTAL variation explained
      model.var(j)    = xLabels(idx);  % variable label
      
      % Prepare for next iteration:
      W             = [W X(:,g(idx,:))]; % selected variable becomes covariable
      X(:,g(idx,:)) = [];            % remove covariable from pool
      uGrp(idx)     = [];            % variable group removed from pool
      delCol        = logical(g(idx,:)==1); % find associated columns to remove
      g(:,delCol)   = [];            % remove associated column
      g(idx,:)      = [];            % indices to variable group removed from pool
      xLabels(idx)  = [];            % label removed from pool
      clear marg;                    % clean up for next iteration
 
   else % Don't add anymore variables to model
      break
   end
end
% =========================================================================



% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('                Stepwise db-RDA:                 \n');
   fprintf('--------------------------------------------------\n');

   fprintf('Global Test: (all variables included)\n');
   headerCell = {'F','p ','R2','R2adj'};    % set up column labels
   resultCell = f_struct2flat(glob);        % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('\nUsed %d permutation of residuals under a FULL model\n',iter(1));
   fprintf('--------------------------------------------------\n\n');

   fprintf('Conditional Tests: (each variable separately)\n');
   headerCell = {'F','p ','R2','R2adj','Variable'}; % set up column labels
   resultCell = f_struct2flat(cond);              % flatten struct to cell array
   resultCell = [headerCell' resultCell']';       % concat
   disp(resultCell);
   fprintf('\nUsed %d permutation of residuals under a FULL model\n',iter(2));
   fprintf('--------------------------------------------------\n\n');

   fprintf('Marginal Tests: (sequential variable addition)\n');
   headerCell = {'Partial F','p','Partial R2','Partial R2adj','Cum R2adj','Variable'}; % set up column labels
   resultCell = f_struct2flat(model);             % flatten struct to cell array
   resultCell = [headerCell' resultCell']';       % concat
   disp(resultCell);
   fprintf('\nUsed %d permutation of residuals under a REDUCED model\n',iter(1));
   if (stop==1) % Stopping rules enforced:
      fprintf('(Only variables that significantly contribute to the model are shown)\n');
   else
      fprintf('(ALL variables are shown, whether significant or not)\n');
   end
   fprintf('--------------------------------------------------\n\n');
   
   fprintf('        R2    = fraction of total variance explained. \n');
   fprintf('        R2adj = adjusted R2. \n');
   fprintf('   Partial R2 = fraction of variance explained after effects of\n');
   fprintf('                variables already in model have been removed. \n');
   fprintf('Partial R2adj = adjusted Partial R2\n');
   fprintf('    Cum R2adj = cumulative fraction of adjusted total variance explained. \n');
   fprintf('--------------------------------------------------\n\n');

   % Report which stopping criteria were met:
   if (stop==1) && (~aFlag || ~rFlag)
      if (~aFlag && ~rFlag)
         fprintf('Variable addition halted due to: ALPHA LEVEL & R^2adj \n');
      elseif (~aFlag)
         fprintf('Variable addition halted due to: ALPHA LEVEL \n');
      else
         fprintf('Variable addition halted due to: R^2adj \n');
      end
   end
end
% ---------------------------------
