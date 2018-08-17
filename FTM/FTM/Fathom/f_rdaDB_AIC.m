function [best,cond] = f_rdaDB_AIC(yDis,ncY,X,BIC,verb,cut,xLabels)
% - stepwise forward selection of explanatory variables in db-RDA using AIC (or BIC)
%
% USAGE: [best,cond] = f_rdaDB_AIC(yDis,ncY,X,BIC,verb,cut,xLabels);
%
% yDis    = square symmetric distance matrix derived from response variables
% ncY     = # columns of original transformed data used to derive yDis
% X       = matrix of explanatory variables
% BIC     = use Bayesian Information Criterion vs. AIC             (default = 0)
% verb    = optionally send result to display                      (default = 1)
% cut     = cut-off value of delta for 'none' indicating insubstantial evidence
%           exists for further addition of variables               (default = 2)
% xLabels = cell array of X labels; if empty, autocreate
%           (e.g., xLabels = {'temp' 'sal' 'depth'})
%
%
% model = structure of marginal effects of the selected model:
%  .RSS   = residual sum-of-squares
%  .R2    = CUMULATIVE fraction of total variance explained
%  .R2adj = CUMULATIVE adjusted R2
%  .AIC   = corrected AIC (or BIC)
%  .var   = variable labels
%  .idx   = index to the selected variables
%  .null  = value of AIC for a null model (i.e., X is only an intercept term)
%  .X     = matrix of explanatory variables of the selected model
%
% cond = structure of conditional effects with the following fields:
%  .RSS   = residual sum-of-squares
%  .R2    = fraction of total variance explained
%  .R2adj = cumulative adjusted R2
%  .AIC   = corrected AIC (or BIC)
%  .var   = variable labels
%
% SEE ALSO: f_rdaDB, f_rdaAIC, f_AIC, f_rdaStepwise

% -----Notes:-----
% This function is used to perform stepwise selection of explanatory variables
% in db-RDA via forward addition based on AIC (or BIC). AIC is used to estimate
% Kullback-Leibler information loss by having a 'lack-of-fit' term and a '# of
% parameters' penalty. An optimal subset of variables, in terms of parsimony, is
% thus achieved by minimizing AIC (= maximizing Akaike weights). In terms of
% AIC, the 'best' model is the one most supported by the empirical data.  Note
% that stepwise selection procedures are not guaranteed to find the 'best'
% subset of variables, only an optimal one, since the minimum found for AIC may
% only be a local (vs. global) one. However, AIC-based procedures do not suffer
% from problems associated with methods based on F-ratios or p-values, namely
% multiple comparison tests, arbitrary significance levels, and false null
% hypotheses.
%
% Model selection based on BIC may result in fewer variables than those
% based on AIC.

% -----References:-----
% Anderson, D. R., K. P. Burnham, and W. L. Thompson. 2000. Null hypothesis
%  testing: problems, prevalence, and an alternative. J. Wildl. Manage.
%  64(4): 912-923.
% Burnham, K. P. and D. R. Anderson. 2001. Kullback-Leibler information as
%  a basis for strong inference in ecological studies. Wildlife Research 28:
%  111-119.
% Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial-modelling: a
%   comprehensive framework for principal coordinate analysis of neighbor
%   matrices (PCNM). Ecological Modelling 196: 483-493.
% Godinez-Dominguez, E. and J. Freire. 2003. Information-theoretic approach
%  for selection of spatial and temporal models of community organization.
%  Mar. Ecol. Prog. Ser. 253: 17-24.

% -----Author:-----
% by David L. Jones, July-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% modified after f_rdaAIC to support db-RDA
 
% Nov-2013: updated documentation; fixed column labels for output table: 'wts'
%           and 'delta' were in the wrong order.

% -----Set defaults & check input:-----
if (nargin < 4), BIC      =   0; end % return AIC by default
if (nargin < 5), verb     =   1; end % send output to display by default
if (nargin < 6), cut      =   2; end % default cut-off value of 2
if (nargin < 7), xLabels  =  cellstr(num2str([1:size(X,2)]'))'; end % default X labels
   
% If labels are not cell arrays, try forcing them:
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end;
xLabels   = xLabels(:); % force row vector
allLabels = xLabels;    % keep a master copy

% Check input:
if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix!');
end

% Check input sizes:
ncX = size(X,2);    % # explanatory variables
n   = size(yDis,1); % # of observations

if n   ~= size(X,1), error('yDis & X need same # of rows!'); end;
if ncX ~= size(xLabels,1)
   error('# of xLabels and explanatory variables do not match!');
end
% -------------------------------------

% Add label specifying no variable addition:
xLabels(ncX+1) = {'none'};

% Preallocate:
cond.RSS(ncX+1,1)    = NaN;
cond.R2(ncX+1,1)     = NaN;
cond.R2adj(ncX+1,1)  = NaN;
cond.AIC(ncX+1,1)    = NaN;
cond.wt(ncX+1,1)     = NaN;
cond.delta(ncX+1,1)  = NaN;
cond.ratio(ncX+1,1)  = NaN;
cond.var(ncX+1)    = {'NaN'};
%
model.RSS(ncX+1,1)   = NaN;
model.R2(ncX+1,1)    = NaN;
model.R2adj(ncX+1,1) = NaN;
model.AIC(ncX+1,1)   = NaN;
model.wt(ncX+1,1)    = NaN;
model.delta(ncX+1,1) = NaN;
model.var(ncX+1)   = {'NaN'};
%
deltaN(ncX+1,1)      = NaN;

% =========================================================================
% CONDITIONAL EFFECTS:
% fprintf('\nPerforming CONDITIONAL TESTS on %d variables...\n',ncX);

k = 1; % 1 variable at a time for conditional tests
for i = 1:ncX
   result        = f_rdaDB(yDis,ncY,X(:,i),[0],0,0,1,0);
   cond.RSS(i)   = result.SSe; % Residual Sum-of-Squares
   cond.R2(i)    = result.R2;
   cond.R2adj(i) = result.R2adj;
   cond.AIC(i)   = f_AIC(cond.RSS(i),n,k+1,BIC); % AIC (corrected)
   cond.var(i)   = xLabels(i);                   % variable label
   
   % Get AIC for a null model (i.e., only an intercept):
   if (i==1)
      cond.RSS(ncX+1)   = result.SSt; % TOTAL SS
      cond.R2(ncX+1)    = NaN;
      cond.R2adj(ncX+1) = NaN;
      cond.AIC(ncX+1)   = f_AIC(cond.RSS(ncX+1),n,k,BIC); % AIC (corrected)
      cond.var(ncX+1)   = {'none'};                       % variable label
   end
   
   % Clean up:
   clear result;
end

% Akaike weights:
cond.delta = cond.AIC - min(cond.AIC); % AIC differences between models
cond.wt    = exp(-0.5*cond.delta);     % relative likelihood of each model
cond.wt    = cond.wt./sum(cond.wt);    % normalize, so sum to 1

% Determine the single best variable:
idx = find(cond.wt == max(cond.wt));   % largest wt wins
if (length(idx)>1), idx = idx(1); end; % in case of ties, select one

% Evidence ratio:
cond.ratio = [cond.wt(:)/cond.wt(idx)].^(-1);

% Sort by wt descending:
[nil,key]   = sort(cond.wt);
key         = flipud(key(:));
sCond.RSS   = cond.RSS(key);
sCond.R2    = cond.R2(key);
sCond.R2adj = cond.R2adj(key);
sCond.AIC   = cond.AIC(key);
sCond.delta = cond.delta(key);
sCond.wt    = cond.wt(key);
sCond.ratio = cond.ratio(key);
sCond.var   = cond.var(key);

% Record delta associated with no variable addition:
deltaN(1) = cond.delta(ncX+1);

% Stop if no variable performs better than 'none':
if (idx == ncX+1)
   best = NaN;
   fprintf('\nMarginal tests skipped: no variable is better than a null model!\n\n')
   return
end

% =========================================================================
% MARGINAL EFFECTS:
if (verb>0)
   fprintf('\nPerforming MARGINAL tests on %d variables...\n',ncX);
end

% Initialize model with best variable from conditional test:
model.RSS(1)   = cond.RSS(idx);
model.R2(1)    = cond.R2(idx);
model.R2adj(1) = cond.R2adj(idx);
model.AIC(1)   = cond.AIC(idx);
model.wt(1)    = cond.wt(idx);
model.delta(1) = cond.delta(idx);
model.var(1)   = cond.var(idx);

W            = X(:,idx); % add best variable to W
X(:,idx)     = [];       % remove best variable from pool
xLabels(idx) = [];       % label removed from pool

clear idx;

for j = 2:ncX            % 1st variable selected from conditional tests
   % fprintf('  -> selecting variable %d of %d \n',j,ncX);
   % Marginal Tests:
   nc_x = size(X,2);     % update # of remaining variables
   k    = size(W,2) + 1; % # variables in [W X(:,i)]
   
   for i = 1:nc_x
      result        = f_rdaDB(yDis,ncY,[W X(:,i)],[0],0,0,1,0);
      marg.RSS(i)   = result.SSe;                   % Residual Sum-of-Squares
      marg.R2(i)    = result.R2;
      marg.R2adj(i) = result.R2adj;
      marg.AIC(i)   = f_AIC(marg.RSS(i),n,k+1,BIC); % AIC (corrected)
      marg.var(i)   = {'none'};                     % variable label
   
      % Get AIC for a model with NO variables added
      if (i==1)
         marg.RSS(nc_x+1)   = model.RSS(j-1); % TOTAL SS
         marg.R2(nc_x+1)    = NaN;
         marg.R2adj(nc_x+1) = NaN;
         marg.AIC(nc_x+1)   = model.AIC(j-1);
         marg.var(nc_x+1)   = {'none'};       % variable label
      end
      
      % Clean up:
      clear result;
   end
   
   % Akaike weights:
   marg.delta = marg.AIC - min(marg.AIC); % AIC differences between models
   marg.wt    = exp(-0.5*marg.delta);     % relative likelihood of each model
   marg.wt    = marg.wt./sum(marg.wt);    % normalize, so sum to 1
    
   % Determine the single best variable:
   idx = find(marg.wt == max(marg.wt));   % largest wt wins
   if (length(idx)>1), idx = idx(1); end  % in case of ties, select one

   % Record delta associated with no variable addition:
   deltaN(j) = marg.delta(nc_x+1);
   
   % Determine if there's sufficient evidence for variable addition:
   if (deltaN(j) <= cut), idx = nc_x+1; end % 'none'
   
   % Add this variable to model:
   if (isempty(idx))
      model.RSS(j)   = NaN;    % Residual Sum-of-Squares
      model.R2(j)    = NaN;
      model.R2adj(j) = NaN;
      model.AIC(j)   = NaN;    % AIC (corrected)
      model.delta(j) = NaN;
      model.wt(j)    = NaN;
      model.var(j)   = {'NaN'}; % variable label
   else
      model.RSS(j)   = marg.RSS(idx);   % Residual Sum-of-Squares
      model.R2(j)    = marg.R2(idx);
      model.R2adj(j) = marg.R2adj(idx);
      model.AIC(j)   = marg.AIC(idx);   % AIC (corrected)
      model.delta(j) = marg.delta(idx);
      model.wt(j)    = marg.wt(idx);
      model.var(j)   = xLabels(idx);    % variable label
   end
   
   if (idx==nc_x+1) % stop when adding no variable is the best
      break;
   end
   
   % Prepare for next iteration:
   W            = [W X(:,idx)]; % add best variable to X
   X(:,idx)     = [];           % remove best variable from pool
   xLabels(idx) = [];           % label removed from pool
   clear marg;                  % clean up for next iteration
      
end
% =========================================================================

% Best model
best.RSS    = model.RSS(1:j);
best.R2     = model.R2(1:j);
best.R2adj  = model.R2adj(1:j);
best.AIC    = model.AIC(1:j);
best.wt     = model.wt(1:j);
best.deltaN = deltaN(1:j);
best.var    = model.var(1:j);

% Get idex to variables comprising the best model:
best.idx = [];
 for i=1:j
    best.idx = [best.idx strmatch(best.var{i}, allLabels, 'exact')];
 end

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   if (BIC<1)
      fprintf('AIC-based stepwise forward selection (db-RDA)      \n');
   else
      fprintf('BIC-based stepwise forward selection (db-RDA)      \n');
   end
   fprintf('--------------------------------------------------\n');
   
   fprintf('Conditional Tests: (each variable separately)\n');
   headerCell = {'RSS','R2','R2adj','AIC','delta','wts','ratio','var'}; % column labels
   resultCell = f_struct2flat(sCond);       % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('--------------------------------------------------\n\n');
   
   fprintf('Marginal Tests: (sequential variable addition)\n');
   headerCell = {'RSS','R2','R2adj','AIC','wts','deltaN','var','idx'}; % column labels
   resultCell = f_struct2flat(best);        % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('--------------------------------------------------\n\n');
   fprintf('RSS    = residual sum-of-squares \n');
   fprintf('R2     = fraction of total variance explained \n');
   fprintf('R2adj  = fraction of adjusted total variance explained \n');
   if (BIC==0)
      fprintf('AIC    = corrected AIC \n');
   else
      fprintf('BIC    = corrected BIC \n');
   end
   fprintf('deltaN = delta associated with NO variable addition \n');
   fprintf('wts    = AIC weights \n');
   fprintf('var    = variable labels \n');
   fprintf('idx    = index to selected variables \n');
   fprintf('\n(Note: RSS, R2, and R2adj in Marginal tests are CUMULATIVE) \n\n');
end


% -----Append fields AFTER outputting to display:-----
% model.maxWt = maxWt;
% model.X      = W(:,1:idxMaxWt);
