function result = f_capMLE(B,dis,grp,txt,M,sm,m,verb,iter,tol)
% - maximum likelihood estimator of mixture proportions based on CAP
% 
% USAGE: result = f_capMLE(B,dis,grp,txt,M,sm,m,verb,iter,tol);
%
% B    = baseline data (= TRAINING SET)           (rows = obs, cols = variables)
% dis = dissimilarity measure to apply to B
%       (e.g., dis = 'bc'; see help for f_dis)
% grp  = column vector of integers specifying group membership for rows of B
% txt  = cell array of group labels                    (if empty, autocreate)
%        e.g., txt = {'North' 'South' 'West'};
% M    = mixture data  (= UNKNOWN SET)            
% sm   = use spatial median instead of centroid
% m    = # PCoA axes to retain                                (see f_capOptimal)
% verb = optionally send result to display                         (default = 1)
% iter = max number of iterations for EM algorithm               (default = 100)
% tol  = convergence tolerance                                 (default = 10^-5)
%  
% 
% result = structure of results with the following fields:
%  .theta = proportion of M that are members of each of B's groups
%  .grpM  = integer specifying group membership predicted by MLE
%  .PP    = corresponding posterior probabilities
%  .nIter = # iterations required to reach convergence
%  .raw   = raw estimates of mixing proportions
%  .PPcap = posterior probabilities from CAP model
%  .marg  = relative measure of certainty of CAP predictions
%  .grpW =  group membership predicted by CAP
% 
% SEE ALSO: f_cap, f_capOptimal, f_mle

% -----Author:-----
% by David L. Jones, May-2013
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 4), txt  = [];    end % default autocreate variable labels
if (nargin < 8), verb = 1;     end % default send results to display
if (nargin < 9), iter = 100;   end % default 100 EM iterations
if (nargin < 10), tol  = 10^-5; end % default tolerance for convergence

grp = grp(:); % force as column vector
if (size(B,1) ~= size(grp,1)), error('Size mismatch b/n B & GRP'); end
if (size(B,2) ~= size(M,2)),   error('B & M must have same # of columns!'); end

% Group labels:
if isempty(txt)
   txt = cellstr(num2str(grp(:)));
   noLabels = 1;
else
   noLabels = 0;
   if (~iscell(txt)), error('TXT must be a cell array of text!'); end
end
txt = txt(:); % force as row vector
if (size(grp,1) ~= size(txt,1))
   error('Size mismatch b/n GRP and TXT!');
end
% -------------------------------------

uGrp  = f_unique(grp);   % unique groups, unsorted
uTxt  = f_unique(txt );  % unique group labels, unsorted
K     = numel(uGrp);     % # of groups
nr    = size(M,1);       % # of unknown observations
theta = repmat(1/K,1,K); % non-informative PRIORS

% Classify M using CAP constructed from B:
cap = f_cap(B,dis,grp,M,sm,0,0,m,0);

% Probability each obs in M is from each of B's group (CAP's distance-based PP):
XIrJ = cap.PP;

%-----Raw estimates:-----
% # obs in M assigned to each of B's groups:
Yi = sum(XIrJ - repmat(max(XIrJ,[],2),1,size(XIrJ,2))==0);

% Raw estimates of proportions:
raw = Yi / sum(Yi);
% -----------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Direct Maximum Likelihood Estimator:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Likelihood each obs is from each group (using this theta):
L = XIrJ .* repmat(theta,nr,1);

% Calculate objective function (= log-likelihood of entire M):
logLik = sum( log(sum(L,2)) ); % Millar, 1987 eq. 7

% -----Expectation-Maximization (EM) Algorithm:-----
dif    = 1; % initialize
nIter  = 0;

while(dif>tol)
   % Keep previous logLik:
   logLikP = logLik;
   
   % Posterior Probability:
   PP = L ./ repmat(sum(L,2),1,K);
   
   % Update theta:
   theta = mean(PP); % (eq 7 in Millar, 1987)
   
   % Update log-likelihood:
   L      = XIrJ .* repmat(theta,nr,1);
   logLik = sum( log(sum(L,2)) );
   
   % Monitor convergence of objective function:
   dif = (logLikP - logLik)/logLikP;
   
   % Monitor # iterations:
   nIter = nIter + 1;
   if (nIter == iter)
      fprintf('\n EM algorithm hasn''t converged after %d iterations!\n',iter)
      result = NaN;
      return
   end
end
% --------------------------------------------------

% Classify unknowns:
[null,grpM] = max(PP,[],2);
% Make sure groups match original labels:
for i = 1:K
   idx = (grpM==uGrp(i));
   grpM(idx) = uGrp(i);
end

% Wrap results up into a structure:
result.theta = theta(:)'; % MLE mixing proportions
result.grpM  = grpM;      % integer specifying predicted group membership
result.PP    = PP;        % corresponding Posterior Probabilities
result.grp   = uGrp(:)';  % list of groups
result.txt   = uTxt(:)';  % cell array of corresponding text labels
result.nIter = nIter;     % # iterations required for convergence
result.raw   = raw;       % raw estimates of mixing proportions
result.PPcap = cap.PP;    % posterior probabilities from CAP model
result.marg  = cap.marg;  % relative measure of certainty of CAP predictions
result.grpW  = cap.grpW;  % group membership predicted by CAP

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf(  '     CAP - MAXIMUM LIKELIHOOD ESTIMATION\n');
   fprintf(  '         Mixing Proportions:\n');
   fprintf(  '--------------------------------------------------\n');
   for i = 1:K
      if (noLabels>0)
         fprintf('Group %d = %-3.4f \n', result.grp(i), result.theta(i));
      else
         fprintf('Group %d = %-3.4f :%s \n', result.grp(i), result.theta(i), result.txt{i});
      end
   end
   fprintf(  '\n# iterations for EM convergence = %d \n',nIter);
   fprintf(  '==================================================\n\n');
end
% ---------------------------------
