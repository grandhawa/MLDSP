function result = f_rdaDB_manova2(yDis,X,Z,iter,model)
% - utility function called by f_rdaDB_manova
%
% USAGE: result = f_rdaDB_manova2(yDis,X,Z,iter,model);
%
% yDis  = square symmetric distance matrix derived from response variables
% X     = col vector of integers specifying levels of Factor 1 for objects in yDis
% Z     = col vector of integers specifying levels of Factor 2 for objects in yDis
% iter  = # iterations for permutation test                        (default = 0)
% model = specifies factor type (fixed = 1, random = 0)

% -----Author:-----
% by David L. Jones, Oct-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Create ANOVA design matrices:
Rx_X         = f_xMatrix(X,1);           % Factor 1
Rx_Z         = f_xMatrix(Z,1);           % Factor 2
[null,Rx_XZ] = f_xMatrix([Rx_X Rx_Z],1); % Interaction 1 x 2
clear null;

% Set up:
n   = size(yDis,1); % # sites
I   = eye(n,n);
uno = ones(n,1);

% Response variable:
G   = (I-(1/n)*(uno*uno'))*(-0.5*(yDis.^2))*(I-(1/n)*(uno*uno')); % Gower's centered matrix
SSt = trace(G); % SS total

if (iter<1)
   iter = 1; % do at least once
else
   Fstat.A  = zeros(iter,1); % preallocate result array
   Fstat.B  = zeros(iter,1);
   Fstat.AB = zeros(iter,1);
end

% ==============================================================================
for i = 1:iter % observed value is considered a permutation
   if (i==1)
      % Use observed G for 1st iteration:
      Gvar_A  = G;
      Gvar_B  = G;
      Gvar_AB = G;
   else
      % Permutation of residuals under a reduced model (L&L,1998
      % p.609-610); add permuted residuals back to form new G:
      Gvar_A  = A.Gfit_W  + f_shuffle(A.Gres_W);
      Gvar_B  = B.Gfit_W  + f_shuffle(B.Gres_W);
      Gvar_AB = AB.Gfit_W + f_shuffle(AB.Gres_W);
   end
   
   % Partition the variation:
   Aperm  = sub_dbRDA(Gvar_A,Rx_X,[Rx_Z Rx_XZ]);  % Factor 1
   Bperm  = sub_dbRDA(Gvar_B,Rx_Z,[Rx_X Rx_XZ]);  % Factor 2
   ABperm = sub_dbRDA(Gvar_AB,Rx_XZ,[Rx_X Rx_Z]); % Factor 1 x 2
   
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

% Wrap results up into a structure:
result(1).so = {'factor 1'};
result(2).so = {'factor 2'};
result(3).so = {'factor 1x2'};
result(4).so = {'residual'};
result(5).so = {'total'};

result(1).df = A.dfR;  % df factor 1
result(2).df = B.dfR;  % df factor 2
result(3).df = AB.dfR; % df interaction
result(4).df = A.dfE;  % df error
result(5).df = n-1;    % df total

result(1).SS = A.SSr;  % SS factor 1
result(2).SS = B.SSr;  % SS factor 2
result(3).SS = AB.SSr; % SS interaction
result(4).SS = A.SSe;  % SS error
result(5).SS = SSt;    % SS total

result(1).MS = A.MSr;  % MS factor 1
result(2).MS = B.MSr;  % MS factor 2
result(3).MS = AB.MSr; % MS interaction
result(4).MS = A.MSe;  % MS error
result(5).MS = NaN;

result(1).F  = Fstat.A(1);  % F-ratio factor 1
result(2).F  = Fstat.B(1);  % F-ratio factor 2
result(3).F  = Fstat.AB(1); % F-ratio interaction
result(4).F  = NaN;
result(5).F  = NaN;

result(1).p = A.p;     % p-value factor 1
result(2).p = B.p;     % p-value factor 2
result(3).p = AB.p;    % p-value interaction
result(4).p = NaN;
result(5).p = NaN;

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
