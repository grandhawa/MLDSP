function result = f_mleCook(M,B,grp,txt,verb,pdf,bs)
% - maximum likelihood estimation via Cook's constrained corrected classification
%
% USAGE: result = f_mleCook(M,B,grp,txt,verb,pdf,bs)
%
% M    = mixture data  (= UNKNOWN SET)            (rows = obs, cols = variables)
% B    = baseline data (= TRAINING SET)
% grp  = column vector of integers specifying group membership for rows of B
% txt  = cell array of variable labels                    (if empty, autocreate)
%           e.g., sLabels = {'Sr' 'Mg' 'Ba'};
% verb = optionally send result to display                         (default = 1)
% pdf  = method of generating probability distribution functions
%        1: multivariate normal + equal group covariance           (default = 1)
%        2: multivariate normal + separate group covariance 
%        3: empirical (distribution-free) based on kernel density estimator
% bs   = # bootstrap iterations for estimates of variance          (default = 0)
%
% result = structure of results with the following fields:
%  .theta = proportion of M that are members of each of B's groups
%  .SD    = bootstrapped standard deviation of theta
%  .cor   = bootstrapped correlations among groups
%  .grpM  = integer specifying predicted group membership of M
%  .PP    = corresponding posterior probabilities
%  .pdf   = method used for generating PDF's
%  .raw   = raw classification proportions
%  .phi   = classification matrix
%
% SEE ALSO: f_mle, f_cda

% -----References:-----
% Portions of this function are based on ideas implemented in Millar's HISEA.exe
% (Windows), Campana & Smith's ISMA.scc (Splus), and White's EM.m (Matlab)
% programs.
%
% Campana, S. E., G. Chouinard, J. M. Hanson, and A. Fréchet. 1999. Mixing
%  and migration of overwintering Atlantic cod (Gadus morhua) stocks near the
%  mouth of the Gulf of St. Lawrence. Can. J. Fish. Aquat. Sci. 56:
%  1873-1881.
% Cook, R. C. 1982. Stock identification of sockeye salmon (Oncorhynchus nerka)
%  with scale pattern recognition. Can. J. Fish. Aquat. Sci. 39: 611-617.
% Cook, R. C. 1983. Simulation and application of stock composition estimators.
%  Can. J. Fish. Aquat. Sci. 40: 2113-2118.
% Millar, R. B. 1987. Maximum likelihood estimation of mixed stock fishery
%  composition. Can J Fish Aquat Sci 44:583-590
% Millar, R. B. 1990. Comparison of methods for estimating mixed stock fishery
%   composition.
% White, J. W. and B. I. Ruttenberg. 2007. Discriminant function analysis in
%  marine ecology: some oversights and their solutions. Mar. Ecol. Prog. Ser.
%  329: 301-305.
%
% For lsqnonneg: http://www.physicsforums.com/showthread.php?t=213313

% -----Author:-----
% by David L. Jones, Feb-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and is
% released under the GNU General Public License, version 2.

% Jun-2012: updated documentation

% -----Check input & set defaults:-----
if (nargin < 4), txt  = []; end % default autocreate variable labels
if (nargin < 5), verb = 1;  end % default send results to display
if (nargin < 6), pdf  = 1;  end % default normal + equal covariance
if (nargin < 7), bs   = 0;  end % default no boostrap iterations

grp = grp(:); % force as column  vector
if (size(B,1) ~= size(grp,1)), error('Size mismatch b/n B & GRP'); end
if (size(B,2) ~= size(M,2)),   error('B & M must have same # of columns!'); end

% Variable labels:
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
nrB     = size(B,1);       % # of baseline observations
priors  = repmat(1/K,1,K); % non-informative PRIORS

if (any(pdf==[1 2])) % --PARAMETRIC:--
   % Get Means & Covariance matrices of BASELINE data:
   mu      = nan(K,nc);   % preallocate
   sigma   = nan(nc,nc,K);
   for i = 1:K
      idx = find(grp == uGrp(i));         % get index to rows of this group
      mu(i,:) = mean(B(idx,:));           % group mean
      if (pdf==1)
         sigma(:,:,i) = f_covPool(B,grp); % pooled covariance matrix
      else
         sigma(:,:,i) = cov(B(idx,:));    % group covariance matrix
      end
   end
   
   % Probability each obs in M is from each of B's group:
   XIrJ = nan(nr,K); % preallocate
   for i = 1:K
      XIrJ(:,i) = mvnpdf(M, mu(i,:), sigma(:,:,i)); % multivariate normal PDF
   end
   
else % --NONPARAMETRIC:--
   % Probability each obs in M is from each of B's group:
   XIrJ = nan(nr,nc,K); % preallocate (row = obs, col = variable, page = group)
   for i = 1:K
      idx = (grp == uGrp(i));
      for j = 1:nc
         XIrJ(:,j,i) = f_empPDF(B(idx,j),M(:,j),0); % empirical PDF
      end
   end
   XIrJ = squeeze(prod(XIrJ,2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  METHOD 1: Constrained Linear Regression (Millar, 1987: p.586):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xij = nan(nrB,K); % preallocate

% -----Estimate PHI using Leave-One-Out (Cook, 1982):-----
for omit = 1:nrB
   % fprintf('%d\n', nrB-omit); % monitor LOO iteration
   B_loo         = B;   % make a copy
   grp_loo       = grp;
   B_loo(omit,:) = [];  % omit this observation
   grp_loo(omit) = [];
   
   if (any(pdf==[1 2])) % --PARAMETRIC:--
      % Get Means & Covariance matrices of BASELINE data (excluding omitted obs):
      mu      = nan(K,nc);   % preallocate
      sigma   = nan(nc,nc,K);
      for i = 1:K
         idx = find(grp_loo == uGrp(i));             % get index to rows of this group
         mu(i,:) = mean(B_loo(idx,:));               % group mean
         if (pdf==1)
            sigma(:,:,i) = f_covPool(B_loo,grp_loo); % single pooled covariance matrix
         else
            sigma(:,:,i) = cov(B_loo(idx,:));        % separate group covariance matrix
         end
      end
      
      % Probability omitted obs in B is from each of B's group (= Xij):
      for i = 1:K
         Xij(omit,i) = mvnpdf(B(omit,:), mu(i,:), sigma(:,:,i)); % multivariate normal PDF
      end
      
   else % --NON-PARAMETRIC:--
      % Probability omitted obs in B is from each of B's group (= Xij):
      Xij_loo = nan(1,nc,K); % preallocate (row = obs, col = variable, page = group)
      for i = 1:K
         idx = (grp_loo == uGrp(i));
         for j = 1:nc
            Xij_loo(1,j,i) = f_empPDF(B_loo(idx,j),B(omit,j),0); % empirical PDF
         end
      end
      Xij(omit,:) = squeeze(prod(Xij_loo,2));
   end
end

% Classify obs to group associated with maximum probability:
grpB = sum( (repmat(1:K,nrB,1) .* ...
   (Xij - repmat(max(Xij,[],2),1,size(Xij,2))==0)) ,2);

% Make sure grpB & grp use same labels:
grpB_old = grpB; % save a copy
for i=1:K
   grpB(grpB_old==i) = uGrp(i);
end
clear grpB_old;

% Classificaion matrix (= PHI):
err = f_errRate(grp,grpB);
phi = err.conf; % Millar, 1990 p. 2238
% clear Xij grpB err;
% --------------------------------------------------------

% # obs in M assigned to each group:
Yi = sum(XIrJ - repmat(max(XIrJ,[],2),1,size(XIrJ,2))==0);

% Raw estimates of proportions:
raw = Yi / sum(Yi);

% Millar's constrained estimates (eq. 12 in Millar, 1987):
% con = (nr*phi')\Yi'; % unconstrained, coefficients can be negative
theta = lsqnonneg((nr*phi'),Yi',priors'); % coefficients always non-negative
theta = theta(:)';

% Likelihood each obs in M is from each of B's groups (using this theta):
L = XIrJ .* repmat(theta,nr,1);

% Posterior Probability:
PP = L ./ repmat(sum(L,2),1,K);

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
      boot           = f_mleCook(M(idxB,:),B,grp,txt,0,pdf,0);
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
result.theta = theta;    % mixture proportions via constrained regression
result.SD    = SD;       % bootstrapped standard deviation of theta
result.cor   = cor;      % bootstrapped correlation among groups
result.grpM  = grpM;     % integer specifying predicted group membership
result.PP    = PP;       % corresponding Posterior Probabilities
result.grp   = uGrp(:)'; % list of groups
result.txt   = uTxt(:)'; % cell array of corresponding text labels
result.raw   = raw(:)';  % raw classification proportions
result.phi   = phi;      % classification matrix

% -----Send output to display:-----
if (verb>0)
   fprintf('\n========================================================\n');
   fprintf(  ' COOK''S CONSTRAINED CORRECTED CLASSIFICATION ESTIMATOR \n');
   fprintf(  '            Mixing Proportions:\n');
   fprintf(  '--------------------------------------------------------\n');
   for i = 1:K
      if (noLabels>0)
         fprintf('Group %d = %-3.4f (%-3.4f) \n', result.grp(i), result.theta(i), result.SD(i));
      else
         fprintf('Group %d = %-3.4f (%-3.4f) :%s \n', result.grp(i), result.theta(i), result.SD(i), result.txt{i});
      end
   end
   fprintf(  '\nPDF = %s\n', pdfTxt);
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
      fprintf('the Cook''s procedure considers this pair similar\n')
   end
   fprintf(  '==================================================\n\n');
end
% ---------------------------------
