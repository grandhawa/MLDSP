function result = f_rdaDB_manova1(yDis,X,W,iter)
% - utility function called by f_rdaDB_manova
%
% USAGE: result = f_rdaDB_manova1(yDis,X,W,iter);
%
% yDis = square symmetric distance matrix derived from response variables
% X    = col vector of integers specifying levels of Factor 1 for objects in yDis
% W    = optional matrix of covariables                               (0 = none)
% iter = # iterations for permutation test                         (default = 0)

% -----Author:-----
% by David L. Jones, Oct-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Create ANOVA design matrices:
X = f_xMatrix(X,1);
if ~isequal(W,0)
   W = f_xMatrix(W,1);
end

% Set up:
n   = size(yDis,1); % # sites
I   = eye(n,n);
uno = ones(n,1);

% Check for covariables:
[nrW,ncW] = size(W);
if (nrW==1) && (ncW==1) && (sum(W(:))==0)
   partial = 0; % no partial analysis
   ncW     = 0; % # covariables
else
   partial = 1; % do partial analysis
   if (ncW==size(X,2))
      if (sum(sum(X==W)) == n)
         error('X & W cannot be the same matrix!')
      else
         W = f_stnd(W); % standardize covariables
      end
   end
end
% -------------------------------------

% Response variable:
G   = (I-(1/n)*(uno*uno'))*(-0.5*(yDis.^2))*(I-(1/n)*(uno*uno')); % Gower's centered matrix
SSt = trace(G); % SS total

if (iter<1)
   iter = 1; % do at least once
else
   Fstat = zeros(iter,1); % preallocate result array
end

% ==============================================================================
for i = 1:iter % observed value is considered a permutation
   if (i==1)
      % Use observed G for 1st iteration:
      Gvar = G;
   else
      % Permutation of residuals under a reduced model (L&L,1998
      % p.609-610); add permuted residuals back to form new G:
      Gvar = Gfit_W  + f_shuffle(Gres_W);
   end
   
   % Partition the variation:
   Aperm  = sub_dbRDA(Gvar_A,Rx_X,[Rx_Z Rx_XZ]);  % Factor 1
   
   % Save observed values:
   if (i==1)
      A  = Aperm;
      B  = Bperm;
      AB = ABperm;
   end
      
   % F-ratios:
   switch model
      case 21 % Factors 1 & 2 are fixed:
         Fstat.A(i)  = Aperm.MSr/Aperm.MSe;
         Fstat.B(i)  = Bperm.MSr/Bperm.MSe;
         Fstat.AB(i) = ABperm.MSr/ABperm.MSe;         
      case 22 % Factor  1 is fixed or random, Factor 2 is random
         Fstat.A(i)  = Aperm.MSr/ABperm.MSr;
         Fstat.B(i)  = Bperm.MSr/ABperm.MSr;
         Fstat.AB(i) = ABperm.MSr/ABperm.MSe;         
      otherwise
         error('Unsupported ANOVA model');
   end
end

% Get permutation-based p-values:
if (iter==1)
   A.p  = NaN;
   B.p  = NaN;
   AB.p = NaN;
else
   % Get permuted stats >= to observed statistic:
   j.A  = find(Fstat.A(2:end)  >=  Fstat.A(1)); % 1st value is observed F
   j.B  = find(Fstat.B(2:end)  >=  Fstat.B(1));
   j.AB = find(Fstat.AB(2:end) >=  Fstat.AB(1));
   
   % Count values & convert to probability:
   A.p  = (length(j.A)+1)./(iter);
   B.p  = (length(j.B)+1)./(iter);
   AB.p = (length(j.AB)+1)./(iter);
end
% ==============================================================================





% =========================================================================
% Permutation Test:
if (iter>0)
   Fperm = zeros(iter-1,1); % preallocate result array
   
   if (partial<1) % NO COVARIABLES:
      % Permutation of residuals under a full model (L&L,1998 p.608):
      for i = 1:(iter-1) % observed value is considered a permutation
         Gperm      = f_shuffle(Gres,2); % just use permuted residuals
         resultPerm = f_rdaDB_manova1(Gperm,X,0,0,1);
         Fperm(i)   = resultPerm.F;
      end
      
   else % WITH COVARIABLES:
      % Permutation of residuals under a reduced model (L&L,1998 p.609-610):
      for i = 1:(iter-1) % observed value is considered a permutation
         Gperm      = Gfit_W + f_shuffle(Gres_W,2); % add permuted residuals back to form new G
         resultPerm = f_rdaDB_manova1(Gperm,X,W,0,1);
         Fperm(i)   = resultPerm.F;
      end
   end
   
   j = find(Fperm >= F);      % get permuted stats >= to observed statistic
   p = (length(j)+1)./(iter); % count values & convert to probability
   
else
   p  = NaN;
end
% =========================================================================




































   
   
   
   
   
   
   
   
   


F = (SSr/ncX)/(SSe/(n-ncX-ncW-1)); % L&L,1998 (eq. 11.19); L,2007 (eq. 3)
% -> this is the same as: F = MSr/MSe

% Coefficient of determination:
SSt   = trace(G);                       % SS total
R2    = SSr/SSt;                        % Legendre,2007 (eq. 5)
R2adj = 1 - ((1-R2)*((n-1)/(n-ncX-1))); % Legendre,2007 (eq. 6)

% If this is a permutation run, stop here:
if (perm>0)
   result.F     = F;
   result.p     = NaN;
   result.R2    = R2;
   result.R2adj = R2adj;
   result.SSe   = SSe;
   result.SSt   = SSt;
   return;
end

% More stats:
dfR = ncX;         % Degrees-of-Freedom regression model
dfE = n-ncX-ncW-1; % Degrees-of-Freedom error
MSr = SSr/dfR;     % Mean Squares regression model
MSe = SSe/dfE;     % Mean Squares error



% -----Wrap results up into a structure:-----
result.F          = F;
result.p          = p;
result.R2         = R2;
result.R2adj      = R2adj;

result.dfR        = dfR;
result.dfE        = dfE;
result.dfT        = n-1;
%
result.SSr        = SSr;
if (partial>0)
   result.SSr_W   = SSr_W;
end
result.SSe        = SSe;
result.SSt        = SSt;

result.MSr        = MSr;
result.MSe        = MSe;
% ------------------------------------------

















%    % =========================================================================
% if (partial<1) % No covariables:
%    [Q1,R1] = qr([uno X],0); H = Q1*Q1';      % Hat-matrix
%    SSr     = trace(H*G*H);                   % SS regression
%    SSe     = trace((I-H)*G*(I-H));           % SS error
%    Gfit    = (H*G*H);                        % fitted values
%    Gres    = G - Gfit;                       % residuals
% else
%    % Variables + covariables:
%    [Q1,R1] = qr([uno X W],0); H_XW = Q1*Q1'; % Hat-matrix
%    SSr_XW  = trace(H_XW*G*H_XW);             % SS regression
%    SSe     = trace((I-H_XW)*G*(I-H_XW));     % SS error
%    
%    % Covariables:
%    [Q1,R1] = qr([uno W],0); H_W = Q1*Q1';    % Hat-matrix for
%    SSr_W   = trace(H_W*G*H_W);               % SS regression for W
%    SSe_W   = trace((I-H_W)*G*(I-H_W));       % SS error for W
%    Gfit_W  = (H_W*G*H_W);                    % fitted values for W
%    Gres_W  = G - Gfit_W;                     % residuals for W
%    % Gres_W = (I-H_W)*G*(I-H_W);             % alternative formulation
%    SSr     = SSr_XW - SSr_W;                 % partial out covariables
% end
% % =========================================================================






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              SUBFUNCTIONS:                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sub = sub_dbRDA(G,X,W)
% - subfunction to partition the variation
%
% USAGE: sub = sub_dbRDA(G,X,W);
% G = Gower's centered matrix
% X = ANOVA design matrix of variable
% W = ANOVA design matrix of covariable

% Setup:
n   = size(G,1);     % # observations
I   = eye(n,n);
uno = ones(n,1);
ncX = size(X,2);     % # columns of X
ncW = size(W,2);     % # columns of W
dfR = ncX;           % df regression model
dfE = (n-ncX-ncW-1); % df error

% =========================================================================
% Variables + covariables:
[Q1,R1] = qr([uno X W],0); H_XW = Q1*Q1'; % Hat-matrix for XW
SSr_XW  = trace(H_XW*G*H_XW);             % SS regression for XW
SSe     = trace((I-H_XW)*G*(I-H_XW));     % SS error for XW

% Covariables:
[Q1,R1] = qr([uno W],0); H_W = Q1*Q1';    % Hat-matrix for W
SSr_W   = trace(H_W*G*H_W);               % SS regression for W
% SSe_W   = trace((I-H_W)*G*(I-H_W));     % SS error for W
Gfit_W  = (H_W*G*H_W);                    % fitted values for W
Gres_W  = G - Gfit_W;                     % residuals for W
SSr     = SSr_XW - SSr_W;                 % partial out covariables
% =========================================================================

% ANOVA stats:
MSr = SSr/dfR; % MS regression
MSe = SSe/dfE; % MS error;

% Wrap result up into a structure
sub.dfR    = dfR;    % df regression
sub.dfE    = dfE;    % df error
sub.SSr    = SSr;    % SS regression
sub.SSe    = SSe;    % SS residuals
sub.MSr    = MSr;    % MS regression
sub.MSe    = MSe;    % MS error
sub.Gfit_W = Gfit_W; % covariable's fitted values
sub.Gres_W = Gres_W; % covariable's residuals




















