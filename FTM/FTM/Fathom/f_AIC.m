function AIC = f_AIC(RSS,n,k,BIC,err)
% - compute AIC (or BIC) for least-squares based methods or classifiers
%
% USAGE: AIC = f_AIC(RSS,n,k,BIC,err);
% 
% RSS = Residual Sum-of-Squares (or classification error rate)
% n   = # observations
% k   = # explanatory variables
% BIC = return BIC instead of AIC                                  (default = 0)
% err = flag indicating RSS specifies a classification error rate  (default = 0)
%                                         
%
% AIC =  AIC (or BIC, corrected for n/k)
%
% SEE ALSO: f_rdaAIC, f_RFaic

% ----- Notes: -----
% This function is used to compute Akaike's Information Criterion (AIC) for
% least-squares based methods for which the residual sum-of-squares is
% provided. When BIC = 1, Schwarz's Bayesian Information Criterion is
% returned instead.

% -----References:-----
% Anderson, D. R., K. P. Burnham, and W. L. Thompson. 2000. Null hypothesis
%  testing: problems, prevalence, and an alternative. J. Wildl. Manage.
%  64(4): 912-923.
% Burnham, K. P. and D. R. Anderson. 2001. Kullback-Leibler information as
%  a basis for strong inference in ecological studies. Wildlife Research 28:
%  111-119.

% -----Author:-----
% by David L. Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2003: k is now NOT adjusted internally to account for intercept
%           and/or error terms
% Feb-2008: updated documentation
% May-2008: added support for out-of-bounds RSS, commented out additive constant
% Sep-2009: added support for classification error rates

%-----Set defaults:-----
if (nargin < 4), BIC = 0; end % return AIC by default
if (nargin < 5), err = 0; end % default RSS specifies RSS


tol = 0.00001; % set tolerance for detecting difference from zero

if (RSS<=tol) % (from 'SpacemakeR for R')
   AIC = NaN;
   return
end

if (err<1)
   % AIC for least-squares based methods:
   if (BIC<1)
      AIC = [n*log(RSS/n)] + (2*k);      % (Anderson et al., 2000)  <- AIC
   else
      AIC = [n*log(RSS/n)] + (log(n)*k); % (from 'R')               <- BIC
   end
else
   % AIC for classification methods (or MSE):
   if (BIC<1)
      AIC = [n*log(RSS)] + (2*k);      %                            <- AIC
   else
      AIC = [n*log(RSS)] + (log(n)*k); %                            <- BIC
   end
end

% Add arbitrary additive constant:
% AIC = [n + n*log(2*pi)] + AIC;       % (from 'R')

% Correct AIC to avoid bias between k and n. This is needed when (n/k)<40,
% but as sample size increases, AIC = AIC_c; so just always use corrected
% version:
AIC = AIC + [(2*k*(k+1)) / (n-k-1)];   % (Anderson et al., 2000)


