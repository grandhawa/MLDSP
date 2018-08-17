function result = f_ancova(Y,grps,W,iter,verb,tol,PW)
% - nonparametric (permutation-based) 1-way ANCOVA
%
% USAGE: result = f_ancova(Y,grps,W,iter,verb,tol)
%
% Y    = column vector of response variable
% grps = column vector of whole numbers specifying group membership        (= Rx)
% W    = column vector specifying covariate
% iter = # iterations for permutation test                          (default = 0)
% verb = optionally send results to display                         (default = 1)
% tol  = minimum p-value needed for covariate effect be significant (default = 0)
%
% structure of results:
%  result.F_slope     = F-stat for test of Rx*Covariate interaction (equal slopes)
%  result.p_slope     = probability that slopes are homogeneous
%  result.F_Rx        = F-stat for test of Rx effect
%  result.p_Rx        = probability of no Rx effect
%  result.F_covariate = F-stat for regression with covariate
%  result.p_covariate = probability of no relationship with covariate
%  result.b           = regression coefficient for covariate
%  result.Ycor        = Y corrected for W
%
% SEE ALSO: f_ancovaPW, aoctool

% -----Notes:-----
% This function is used to perform a nonparametric (permutation-based) Analysis
% of Covariance (ANCOVA).
%
% If 'p_covariate' is at least as small as 'tol', then the response
% variable 'Y' is corrected for the effect of the covariable 'W' and output
% as 'Ycor'. This can also be done outside the function: For cases where
% there is a significant relationship of Y with W, you can correct (remove)
% the effect of W on Y using: Ycorrected = Y - bW where b is the regression
% coefficient for W.

% -----References:-----
% Petraitis, P. S., S. J. Beaupre, & A. E. Dunham. 2001. ANCOVA: Nonparametric
%   and randomization approaches. Pages 116-133 in S. M. Scheiner & J. Gurevitch
%   (eds.), Design and Analysis of Ecological Experiments. Oxford University
%   Press.

% -----Author:-----
% by David L. Jones, Jun-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2010: now initializes RAND to a different state
% Jul-2012: added verbose output

% -----Set defaults & check input:-----
if (nargin < 4), iter = 0; end % no permutation test by default
if (nargin < 5), verb = 1; end % send output to display by default
if (nargin < 6), tol  = 0; end % don't correct for covariable by default
if (nargin < 7), PW   = 0; end % internal flag indicating this is a PW slope run

% Force column vectors:
Y    = Y(:);
grps = grps(:);
W    = W(:);

noRows = size(Y,1);

% Homogeneity of Slopes:
[null,anova,coefs,stats] = aoctool(W,Y,grps,0.5,'','','','off','separate lines');
F_slope = anova{4,5}; % observed F-stat for Homogeneity of Slopes

% Rx Effect and Covariate Regression:
[null,anova,coefs,stats] = aoctool(W,Y,grps,0.5,'','','','off','parallel lines');
F_Rx        = anova{2,5}; % observed F-stat for Rx
F_covariate = anova{3,5}; % observed F-stat for covariate
b           = coefs{size(coefs,1),2}; % collect regression coefficients
% -----------------------------------------------------------------
if (iter>0)
   rand('twister',sum(100*clock)); % initialize RAND to a different state
   
   if PW<0
      fprintf('\nPermuting data %d times...\n',iter-1);
   end
   
   % -----SLOPE:-----
   F_slopePerm     = zeros(iter-1,1); % preallocate result array
   for i = 1:(iter-1) % observed value is considered a permutation
      
      % Permute covariate (Petraitis et al., 2001):
      idx = f_shuffle(1:noRows);
      [null,anovaSlopePerm] = aoctool(W(idx,:),Y,grps,0.5,'','','','off','separate lines');
      F_slopePerm(i)        = anovaSlopePerm{4,5};
   end
   
   j       = find(F_slopePerm >= F_slope); % get permuted stats >= to observed statistic
   p_slope = (length(j)+1)./(iter);        % count values & convert to probability
   % ----------------
   
   % -----Rx and Covariate:-----
   if (PW==0) % not a "PW slope run"
      F_RxPerm        = zeros(iter-1,1);
      F_covariatePerm = zeros(iter-1,1);
      for i = 1:(iter-1) % observed value is considered a permutation
         % Permute rows of Rx only:
         idx = f_shuffle(1:noRows);
         [null,anovaRxPerm] = aoctool(W,Y,grps(idx,:),0.5,'','','','off','parallel lines');
         F_RxPerm(i)        = anovaRxPerm{2,5};
         
         % Permute covariate:
         idx = f_shuffle(1:noRows);
         [null,anovaCovariatePerm] = aoctool(W(idx,:),Y,grps,0.5,'','','','off','parallel lines');
         F_covariatePerm(i)        = anovaCovariatePerm{3,5};
         
      end
      j    = find(F_RxPerm >= F_Rx); % get permuted stats >= to observed statistic
      p_Rx = (length(j)+1)./(iter);  % count values & convert to probability
      
      j            = find(F_covariatePerm >= F_covariate);
      p_covariate  = (length(j)+1)./(iter);
   else % skip because a "PW slope run"
      p_Rx        = NaN;
      p_covariate = NaN;
   end
   % ---------------------------
else
   p_slope     = NaN;
   p_Rx        = NaN;
   p_covariate = NaN;
end
%  -----------------------------------------------------------------

% -----Remove the effect of the covariable:-----
if (p_covariate <= tol)
   Ycor = Y - W*b;
else
   Ycor = NaN;
end
% ----------------------------------------------


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('    Nonparametric (Permutation-based) ANCOVA:\n');
   fprintf('--------------------------------------------------\n');
   fprintf('Rx*Covariate: F = %8.4f  p =  %3.5f \n',F_slope,p_slope);
   fprintf('Rx Effect:    F = %8.4f  p =  %3.5f \n',F_Rx,p_Rx);
   fprintf('Covariate:    F = %8.4f  p =  %3.5f \n',F_covariate,p_covariate);
   fprintf('No. of permutations = %d \n',iter);
   fprintf('==================================================\n');
   fprintf('NOTE:\n');
   fprintf('A significant ''Rx*Covariate interaction'' term suggests\n')
   fprintf('the slopes of the regression lines are NOT homogeneous\n')
end
% ---------------------------------

% Wrap results up into a structure
result.F_slope     = F_slope;
result.p_slope     = p_slope;
result.F_Rx        = F_Rx;
result.p_Rx        = p_Rx;
result.F_covariate = F_covariate;
result.p_covariate = p_covariate;
result.b           = b;
result.Ycor        = Ycor;
