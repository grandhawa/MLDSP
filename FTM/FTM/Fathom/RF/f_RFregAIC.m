function [best,cond] = f_RFregAIC(X,Y,nTree,BIC,verb,cut,xLabels,opt)
% - AIC-based variable selection for Random Forest Regression
%
% USAGE: [best,cond] = f_RFregAIC(X,Y,nTree,BIC,verb,cut,xLabels);
%
% X       = input data                            (rows = obs, cols = variables)
% Y       = column vector of integers specifying class membership
% nTree   = # random trees to grow           (default = 1000 if set to empty [])
% BIC     = use Bayesian Information Criterion                     (default = 0)
% verb    = optionally send result to display                      (default = 1)
% cut     = cut-off value of delta for 'none' indicating insubstantial evidence
%           exists for further addition of variables      (default = 2)
% xLabels = cell array of X labels; if empty, autocreate
%           (e.g., xLabels = {'temp' 'sal' 'depth'})
%
% SEE ALSO: f_RFreg, f_RFaic, f_RFimp, f_AIC, f_rdaAIC

% -----Author:-----
% by David L. Jones, Oct-2011
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), nTree    =  []; end  % use default # trees
if (nargin < 4), BIC      =   0; end; % return AIC by default
if (nargin < 5), verb     =   1; end; % send output to display by default
if (nargin < 6), cut      =   2; end; % default cut-off value of 2
if (nargin < 7), xLabels  =  cellstr(num2str([1:size(X,2)]'))'; end % default X labels
if (nargin < 8), opt.null =  NaN; end % create dummy structure of options

% If labels are not cell arrays, try forcing them:
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end;
xLabels   = xLabels(:); % force row vector
allLabels = xLabels;    % keep a master copy

% Check input sizes:
ncX = size(X,2); % # explanatory variables
n   = size(Y,1); % # of observations

if n   ~= size(X,1), error('Y & X need same # of rows'); end;
if ncX ~= size(xLabels,1)
   error('# of xLabels and explanatory variables do not match!');
end

% Replace empty inputs with defaults:
if isempty(nTree), nTree = 1000; end

% Default no standardization of predictors:
stnd = 'stnd';
% stnd = 'raw';
% -------------------------------------

% Add label specifying no variable addition:
xLabels(ncX+1) = {'none'};

% Preallocate:
cond.RSS(ncX+1,1)    = NaN;
cond.R2(ncX+1,1)     = NaN;
cond.AIC(ncX+1,1)    = NaN;
cond.wt(ncX+1,1)     = NaN;
cond.delta(ncX+1,1)  = NaN;
cond.ratio(ncX+1,1)  = NaN;
cond.var(ncX+1)      = {'NaN'};
%
model.RSS(ncX+1,1)   = NaN;
model.R2(ncX+1,1)    = NaN;
model.AIC(ncX+1,1)   = NaN;
model.wt(ncX+1,1)    = NaN;
model.delta(ncX+1,1) = NaN;
model.var(ncX+1)     = {'NaN'};
%
deltaN(ncX+1,1)      = NaN;

% =========================================================================
% CONDITIONAL EFFECTS:
fprintf('\nPerforming CONDITIONAL TESTS on %d variables...\n',ncX);

k = 1; % 1 variable at a time for conditional tests
for i = 1:ncX
   result      = f_RFreg(X(:,i),Y,nTree,[],0,0,stnd,0,[],opt);
   cond.RSS(i) = result.mse(end)*n;              % Residual Sum-of-Squares
   cond.R2(i)  = result.rsq(end);                % R-squared
   cond.AIC(i) = f_AIC(cond.RSS(i),n,k+1,BIC);   % AIC (corrected)
   cond.var(i) = xLabels(i);                     % variable label
   
   % Get AIC for a null model:
   if (i==1)
      uno             = ones(size(Y));                    % just an intercept term
      result          = f_RFreg(uno,Y,nTree,[],0,0,stnd,0,[],opt);
      cond.RSS(ncX+1) = result.mse(end)*n;                % Residual Sum-of-Square
      cond.R2(ncX+1)  = NaN;
      cond.AIC(ncX+1) = f_AIC(cond.RSS(ncX+1),n,1,BIC);   % AIC (corrected)
      cond.var(ncX+1) = {'none'};                         % variable label
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
model.AIC(1)   = cond.AIC(idx);
model.wt(1)    = cond.wt(idx);
model.delta(1) = cond.delta(idx);
model.var(1)   = cond.var(idx);

W            = X(:,idx); % add best variable to W
X(:,idx)     = [];       % remove best variable from pool
xLabels(idx) = [];       % label removed from pool

clear idx;

for j = 2:ncX            % 1st variable selected from Conditional tests
   fprintf('  -> selecting variable %d of %d \n',j,ncX);
   % Marginal Tests:
   nc_x = size(X,2);     % update # of remaining variables
   k    = size(W,2) + 1; % # variables in [W X(:,i)]
   
   for i = 1:nc_x
      result      = f_RFreg([W X(:,i)],Y,nTree,[],0,0,stnd,0,[],opt);
      marg.RSS(i) = result.mse(end)*n;              % Residual Sum-of-Squares
      marg.R2(i)  = result.rsq(end);                % R-squared
      marg.AIC(i) = f_AIC(marg.RSS(i),n,k+1,BIC);   % AIC (corrected)
      marg.var(i) = xLabels(i);                     % variable label
      
      % Get AIC for a model with NO variables added
      if (i==1)
         marg.RSS(nc_x+1) = model.RSS(j-1); % Residual Sum-of-Squares
         marg.R2(nc_x+1)  = NaN;
         marg.AIC(nc_x+1) = model.AIC(j-1); % AIC (corrected)
         marg.var(nc_x+1) = {'none'};       % variable label
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
      model.RSS(j)   = NaN;   % Residual Sum-of-Squares
      model.R2(j)    = NaN;   % R-squared
      model.AIC(j)   = NaN;   % AIC (corrected)
      model.wt(j)    = NaN;
      model.delta(j) = NaN;
      model.var(j) = {'NaN'}; % variable label
   else
      model.RSS(j)   = marg.RSS(idx);   % Residual Sum-of-Squares
      model.R2(j)    = marg.R2(idx);    % R-squared
      model.AIC(j)   = marg.AIC(idx);   % AIC (corrected)
      model.wt(j)    = marg.wt(idx);
      model.delta(j) = marg.delta(idx);
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
      fprintf('AIC-based stepwise forward selection (RANDOM FOREST)\n');
   else
      fprintf('BIC-based stepwise forward selection (RANDOM FOREST)\n');
   end
   fprintf('--------------------------------------------------\n');
   
   fprintf('Conditional Tests: (each variable separately)\n');
   headerCell = {'RSS','R2','AIC','wts','delta','ratio','var'}; % column labels
   resultCell = f_struct2flat(sCond);       % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('--------------------------------------------------\n\n');
   
   fprintf('Marginal Tests: (sequential variable addition)\n');
   headerCell = {'RSS','R2','AIC','wts','deltaN','var','idx'}; % column labels
   resultCell = f_struct2flat(best);        % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   disp(resultCell);
   fprintf('--------------------------------------------------\n\n');
   fprintf('# Trees = %d\n',nTree);
   fprintf('RSS     = residual sum of squares \n');
   fprintf('R2      = fraction of total variation explained \n');
   if (BIC==0)
      fprintf('AIC     = corrected AIC \n');
   else
      fprintf('BIC     = corrected BIC \n');
   end
   fprintf('deltaN  = delta associated with NO variable addition \n');
   fprintf('wts     = AIC weights \n');
   fprintf('var     = variable labels \n');
   fprintf('idx     = index to selected variables \n');
   fprintf('\n(Note: RSS, R2 in Marginal Tests are CUMULATIVE) \n\n')
end
