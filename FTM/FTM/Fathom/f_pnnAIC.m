function sol = f_pnnAIC(Y,X,sm,BIC,verb,xLabels)
% - stepwise selection of explanatory variables in PNN using AIC (or BIC)
%
% USAGE: sol = f_pnnAIC(grp,X,sm,BIC,verb,xLabels)
%
% Y        = vector of integers specifying group membership
% X        = matrix of explanatory variables
% sm       = smoothing factor
% BIC      = use Bayesian Information Criterion (default = 0)
% verb     = optionally send result to display  (default = 1)
%
% xLabels  = cell array of X labels; if empty, autocreate
%            (e.g., xLabels = {'temp' 'sal' 'depth'})
%
% sol     =   cell-array of results; each cell contains a structure of
%             values for each step of the variable selection process,
%             sorted by sol.wt, descending
% sol.logLik = log likelihood
% sol.noPar  = # estimable parameters
% sol.AIC    = corrected AIC (or BIC)
% sol.delta  = difference between AIC of model and 'best' model
% sol.wt     = Akaike weights, normalized to sum = 1
% sol.ratio  = evidence ratio
% sol.var    = variable added to (or deleted from) model
%
% SEE ALSO: f_AIC

% -----Author:-----
% by David L. Jones, Jan-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), sm       =   0.1; end; % return AIC by default
if (nargin < 4), BIC      =   0;   end; % return AIC by default
if (nargin < 5), verb     =   1;   end; % send output to display by default 
if (nargin < 6), xLabels  = cellstr(num2str([1:size(X,2)]'))'; end; % default X labels


% If labels are not cell arrays, try forcing them:
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end;
xLabels = xLabels(:); % force row vector

% Check input sizes:
ncX = size(X,2); % # explanatory variables
n   = size(Y,1); % # of observations

if n   ~= size(X,1), error('Y & X need same # of rows'); end;
if ncX ~= size(xLabels,1)
   error('# of xLabels and explanatory variables do not match!');
end


% FORWARD ADDITION:  
fprintf('Performing forward addition of: \n');

W = []; % initialize with no variables

for k = 1:ncX
   fprintf('  -> %d variables \n',k);
   
   % Create alternative models by adding 1 variable:
   nc = size(X,2);
   for i = 1:nc;
      temp     = f_pnn([W X(:,i)],Y,[W X(:,i)],'empirical',sm);
      resid    = (temp.target - temp.Bayes);  % residuals
      
      RSS(i)   = sum(sum(resid.^2));           % Residual Sum-of-Squares
      AIC(i)   = f_AIC(RSS(i),n,k+1,BIC);      % AIC (corrected)
      noPar(i) = k+2;                          % # estimable parameters
      added(i) = xLabels(i);                   % label of variable added
   end
   
   % Include a model with NO variables added:
   if (k>1)
      temp        = f_pnn([W],Y,[W],'empirical',sm);
      resid       = (temp.target - temp.Bayes);  % residuals
      
      RSS(nc+1)   = sum(sum(resid.^2));           % Residual Sum-of-Squares
      % RSS(nc+1) = 1 - sum(result.canVar);
   else
      RSS(nc+1)   = 99999999999999999999999; % just intercept + error
   end   
   
   AIC(nc+1)   = f_AIC(RSS(nc+1),n,k+1,BIC);             
   noPar(nc+1) = (k-1) + 2;
   added(nc+1) = {'none'};  % no variable added
   
   % Akaike weights:
   delta = AIC - min(AIC);  % AIC differences between models    
   wt    = exp(-0.5*delta); % relative likelihood of each model 
   wt    = wt./sum(wt);     % normalize, so sum to 1            
   
   % Select model with largest wt (= smallest AIC)
   idx = find(wt == max(wt));
   
   % Sort by wt descending
   [null,key] = sort(wt);
   key        = flipud(key(:));
   RSS        = RSS(key);
   noPar      = noPar(key);
   AIC        = AIC(key);
   delta      = delta(key);
   wt         = wt(key);
   added      = added(key);
   
   % Save for output:   
   sol{k}.logLik = log(RSS(:)/n);      % maximized log-likelihood
   sol{k}.noPar  = noPar(:);
   sol{k}.AIC    = AIC(:);
   sol{k}.delta  = delta(:);
   sol{k}.wt     = wt(:);
   sol{k}.ratio  = [wt(:)/wt(1)].^(-1); % evidence ratio
   sol{k}.var    = added(:);
   
   
   if idx > size(X,2)
      break
   end
   
   % Prepare for next iteration:
   W            = [W X(:,idx)];          % add variable to W
   X(:,idx)     = [];                    % removed from pool of available variables
   xLabels(idx) = [];                    % remove label from pool
   clear RSS noPar AIC delta wt added;   % clean up
   
end



% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');   
   if (BIC<1)
      fprintf('AIC-based forward selection in PNN:\n');
   else
      fprintf('BIC-based forward selection in REDUNDANCY ANALYSIS:\n');
   end
   fprintf('--------------------------------------------------\n');
   
   % Set up column labels:
   if (BIC<1)
      headerCell = {'logLik' 'k' 'AIC' 'delta' 'wts' 'ratio' 'VAR'};
   else
      headerCell = {'logLik' 'k' 'BIC' 'delta' 'wts' 'ratio' 'VAR'};
   end
   
   % Display results for each step separately:
   for i = 1:k
      fprintf('Models with %d variables added: \n',i);
      resultCell = f_struct2flat(sol{i});        % flatten struct to cell array
      resultCell = [headerCell' resultCell']';   % concat
      disp(resultCell);
      fprintf('--------------------------------------------------\n\n');
   end
   
   % Legend:
   fprintf('logLik  = value of maximized log-likelihood \n');
   fprintf('k       = # estimable parameters (includes error term) \n');
   
   if (BIC<1)
      fprintf('AIC     = Akaike Information Criterion \n');
   else
      fprintf('BIC     = Bayesian Information Criterion \n');
   end
   
   fprintf('delta   = difference between AIC and min(AIC) \n');
   
   if (BIC<1)
      fprintf('wts     = Akaike wts (~ probability of being ''best'') \n');
   else
      fprintf('wts     = Bayesian wts \n');
   end
   
   fprintf('ratio   = evidence ratio \n');
   fprintf('VAR     = variable added to model \n');
   fprintf('--------------------------------------------------\n\n');
end








