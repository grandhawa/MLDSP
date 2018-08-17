function result = f_mle(M,B,grp,txt,verb,pdf,iter,bs,tol)
% - maximum likelihood estimator of mixture proportions via an EM algorithm
%
% USAGE: result = f_mle(M,B,grp,txt,verb,pdf,iter,bs,tol)
%
% M    = mixture data  (= UNKNOWN SET)            (rows = obs, cols = variables)
% B    = baseline data (= TRAINING SET)
% grp  = column vector of integers specifying group membership for rows of B
% txt  = cell array of group labels                       (if empty, autocreate)
%        e.g., txt = {'North' 'South' 'West'};
% verb = optionally send result to display                         (default = 1)
% pdf  = method of generating probability distribution functions
%        1: multivariate normal + equal group covariance           (default = 1)
%        2: multivariate normal + separate group covariance
%        3: empirical (distribution-free) based on kernel density estimator
% iter = max number of iterations for EM algorithm               (default = 100)
% bs   = # iterations for bootstrapped estimates of variance       (default = 0)
% tol  = convergence tolerance                                 (default = 10^-5)
%
% result = structure of results with the following fields:
%  .theta = proportion of M that are members of each of B's groups
%  .SD    = bootstrapped standard deviation of theta
%  .cor   = bootstrapped correlations among groups
%  .grpM  = integer specifying predicted group membership of M
%  .PP    = corresponding posterior probabilities
%  .nIter = # iterations required to reach convergence
%  .pdf   = method used for generating PDF's
%  .raw   = raw estimates of mixing proportions
%
% SEE ALSO: f_mleCook, f_cda, f_covPool

% -----References:-----
% Portions of this function are based on ideas implemented in Millar's HISEA.exe
% (Windows), Campana & Smith's ISMA.scc (Splus), and White's EM.m (Matlab)
% programs. This function had been tested against these and provides similar
% results.
%
% Campana, S. E., G. Chouinard, J. M. Hanson, and A. Fréchet. 1999. Mixing
%  and migration of overwintering Atlantic cod (Gadus morhua) stocks near the
%  mouth of the Gulf of St. Lawrence. Can. J. Fish. Aquat. Sci. 56:
%  1873-1881.
% Millar, R. B. 1987. Maximum likelihood estimation of mixed stock fishery
%  composition. Can J Fish Aquat Sci 44:583-590
% Millar, R. B. 1990. Comparison of methods for estimating mixed stock fishery
%   composition.
% White, J. W. and B. I. Ruttenberg. 2007. Discriminant function analysis in
%  marine ecology: some oversights and their solutions. Mar. Ecol. Prog. Ser.
%  329: 301-305.

% -----Note:-----
% The 'log-likelihood' is used instead of the 'likelihood' in this scenario
% for mathematical convenience, as it lets you turn multiplication into
% addition.

% -----Author:-----
% by David L. Jones, Jan-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2011: moved Cook's constrained regression method to 'f_mleCook.m'; added
%           empirical PDF; record the method used
% May-2013: updated documentation; clarified that input parameter 'txt' is
%           a group (not variable) label.

% -----Check input & set defaults:-----
if (nargin < 4), txt  = [];    end % default autocreate variable labels
if (nargin < 5), verb = 1;     end % default send results to display
if (nargin < 6), pdf  = 1;     end % default normal + equal covariance
if (nargin < 7), iter = 100;   end % default 100 EM iterations
if (nargin < 8), bs   = 0;     end % default no boostrap iterations
if (nargin < 9), tol  = 10^-5; end % default tolerance for convergence

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

% Check PDF's:
if (~any(pdf==[1 2 3]))
   error('PDF must be 1, 2, or 3!')
end

% Record the method used:
switch pdf
   case 1
      pdfTxt = 'MVN + equal covariance';
   case 2
      pdfTxt = 'MVN + unequal covariances';
   case 3
      pdfTxt = 'Distribution-free (empirical)';
end
% -------------------------------------

uGrp    = f_unique(grp);   % unique groups, unsorted
uTxt    = f_unique(txt);   % unique group labels, unsorted
K       = numel(uGrp);     % # of groups
[nr,nc] = size(M);         % # of unknown observations, variables
theta   = repmat(1/K,1,K); % non-informative PRIORS

% -----Probability Density Functions (PDF's):-----
if (any(pdf==[1 2])) % --PARAMETRIC:--
   % Get Means & Covariance matrices of BASELINE data:
   mu      = nan(K,nc);   % preallocate
   sigma   = nan(nc,nc,K);
   for i = 1:K
      idx = find(grp == uGrp(i));         % get index to rows of this group
      mu(i,:) = mean(B(idx,:));           % group mean
      if (pdf==1)
         sigma(:,:,i) = f_covPool(B,grp); % single pooled covariance matrix
      else
         sigma(:,:,i) = cov(B(idx,:));    % separate group covariance matrix
      end
   end
   
   % Probability each obs in M is from each of B's group:
   XIrJ = nan(nr,K); % preallocate
   for i = 1:K
      XIrJ(:,i) = mvnpdf(M, mu(i,:), sigma(:,:,i)); % multivariate normal PDF
   end
else
   XIrJ = nan(nr,nc,K); % preallocate (row = obs, col = variable, page = group)
   for i = 1:K
      idx = (grp == uGrp(i));
      for j = 1:nc
         XIrJ(:,j,i) = f_empPDF(B(idx,j),M(:,j),0); % empirical PDF
      end
   end
   XIrJ = squeeze(prod(XIrJ,2));
end
% ------------------------------------------------

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

% -----Bootstrapped estimates of variance in theta (Millar, 1987):-----
if (bs>0)
   thetaBoot = nan(bs,K); % preallocate
   for i=1:bs
      idxB           = randi(nr,nr,1); % get bootstrapped index to rows of M
      boot           = f_mle(M(idxB,:),B,grp,txt,0,pdf,iter,0,tol);
      thetaBoot(i,:) = boot.theta;
      clear boot;
   end
   SD  = std(thetaBoot);  % standard deviation
   cor = corr(thetaBoot); % correlation matrix
else
   SD  = nan(1,K);
   cor = NaN;
end
% ------------------------------------------------------

% Wrap results up into a structure:
result.theta = theta(:)'; % MLE mixing proportions
result.SD    = SD;        % bootstrapped standard deviation of theta
result.cor   = cor;       % bootstrapped correlation among groups
result.grpM  = grpM;      % integer specifying predicted group membership
result.PP    = PP;        % corresponding Posterior Probabilities
result.grp   = uGrp(:)';  % list of groups
result.txt   = uTxt(:)';  % cell array of corresponding text labels
result.nIter = nIter;     % # iterations required for convergence
result.pdf   = pdfTxt;    % method for generating PDF's
result.raw   = raw;       % raw estimates of mixing proportions

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf(  '        MAXIMUM LIKELIHOOD ESTIMATION\n');
   fprintf(  '            Mixing Proportions:\n');
   fprintf(  '--------------------------------------------------\n');
   for i = 1:K
      if (noLabels>0)
         fprintf('Group %d = %-3.4f (%-3.4f) \n', result.grp(i), result.theta(i), result.SD(i));
      else
         fprintf('Group %d = %-3.4f (%-3.4f) :%s \n', result.grp(i), result.theta(i), result.SD(i), result.txt{i});
      end
   end
   fprintf(  '\nPDF = %s\n', pdfTxt);
   fprintf(  '\n# iterations for EM convergence = %d \n',nIter);
   fprintf(  '(Bootstrapped SD''s are in parentheses)\n\n');
   
   if (bs>0)
      fprintf(  '--------------------------------------------------\n');
      fprintf('Bootstrapped correlations among Groups:\n');
      [null,null2,pairList] = f_subsetDisPW(1-result.cor,[1:K]');
      corVec                = f_unwrap(result.cor,0); % unwrap matrix to a vector
      for i = 1:size(pairList,1)
         fprintf('Group %d vs. %d = %+3.4f \n', result.grp(pairList(i,1)), result.grp(pairList(i,2)), corVec(i));
      end
      fprintf('\n# iterations for Bootstrap = %d \n\n',bs);
      fprintf('NOTE: Strong neg. correlations between groups suggests\n')
      fprintf('the MLE procedure considers this pair similar\n')
   end
   fprintf(  '==================================================\n\n');
end
% ---------------------------------

