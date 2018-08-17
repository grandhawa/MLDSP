function result = f_distlmPartial(n,I,G,x,c,iter);
% - utility function called by f_distlm

% - performs a partial DISTLM when co-variables are specified

% -----Input/Output:-----
% n = # rows/colums in distance matrix
% I = I matrix
% G = Gower's centered matrix
% x = matrix of explanatory variables
% c = matrix of co-variables
% iter = # iterations for permutation test

% -----References:-----
% % Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. (Section 10.5)

% -----Author:-----
% by David L. Jones, Aug-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% add intercept:
xx = [ones(n,1) x];     % explanatory matrix
cc = [ones(n,1) c];     % co-variable matrix
xc = [ones(n,1) [x c]]; % combined matrix

% comput Hat-matrix via QR:
[Q1,R1] = qr(xx,0); xH=  Q1*Q1';
[Q1,R1] = qr(cc,0); cH=  Q1*Q1';
[Q1,R1] = qr(xc,0); xcH= Q1*Q1';

xM  = size(xx,2); % # of parameters in explanatory matrix
cM  = size(cc,2); % # of parameters in co-variable matrix

% Sum-of-Squares:
SSR  = trace(xH*G*xH);   % SS Regression (McArdle & Anderson, 2001)
SSC  = trace(cH*G*cH);   % SS co-variables
SSRC = trace(xcH*G*xcH); % SS Regression + co-variables
SST  = trace(G);         % SS Total

% Degrees of Freedom:
df.R = xM-1;               % df Regression
df.C = cM-1;               % df Co-variables
df.T = n-1;                % df Total
df.E = df.T - df.R - df.C; % df Error, corrected for co-variables      

% Factor out effect of co-variables:
SSR = SSRC - SSC;
SSE = (SST - SSR) - SSC;


% Mean Squares & F-Stat:
MSR = SSR/df.R;
MSE = SSE/df.E;
F   = MSR/MSE;

% Permutation of Residuals under the Reduced Model:
if (iter>0)
   fprintf('\nPermuting residuals (reduced model) %d times...\n',iter-1);
   rand('state',sum(100*clock)); % set random generator to new state
   
   % -----Reduced Model:-----
   R       = (I-cH)*G*(I-cH); % residuals (McArdle & Anderson, 2001)
   G_fit   = cH*G*cH;         % fitted values      
   
   F_perm = zeros(iter-1,1); % preallocate results array
   for i=1:iter
      R_perm = f_shuffle(R,2); % permute residuals, preserving 'trace'
      G_perm = G_fit + R_perm; % add permuted residuals back to G-fit         
      
      % Permuted Sum-of-Squares:
      SSR_perm  = trace(xH*G_perm*xH);   % SS Regression
      SSC_perm  = trace(cH*G_perm*cH);   % SS Co-variables
      SSRC_perm = trace(xcH*G_perm*xcH); % SS Regression + co-variables
            
      % Factor out effect of co-variables:
      SSR_perm = SSRC_perm - SSC_perm;
      SSE_perm = (SST - SSR_perm) - SSC_perm;
      
      % Permuted Mean Squares & F-stat:
      MSR_perm  = SSR_perm/df.R;
      MSE_perm  = SSE_perm/df.E;
      F_perm(i) = MSR_perm/MSE_perm;
      
   end
   j = find(F_perm >= F);     % get permuted F >= to observed F
   p = (length(j)+1)./(iter); % count values & convert to probability      
   
else
   p = NaN;   
end

% Wrap results up into a structure:      
result(1).so = {'Co-variables'};
result(2).so = {'Regression'};
result(3).so = {'Residual'};
result(4).so = {'Total'};      

result(1).df = df.C;
result(2).df = df.R;
result(3).df = df.E;
result(4).df = df.T;

result(1).SS = SSC;
result(2).SS = SSR;
result(3).SS = SSE;
result(4).SS = SST;

result(1).MS = NaN;
result(2).MS = MSR;
result(3).MS = MSE;
result(4).MS = NaN;

result(1).F = NaN;
result(2).F = F;
result(3).F = NaN;
result(4).F = NaN;

result(1).p = NaN;
result(2).p = p;
result(3).p = NaN;
result(4).p = NaN;

result(1).R2 = SSC/SST;
result(2).R2 = SSR/SST;
result(3).R2 = SSE/SST;
result(4).R2 = SST/SST;

