function [df,SS,MS,F,p,result] = f_distlmUtil(n,m,H,I,G,iter)
% - utility function called by f_distlm

% -----Input/Output:-----
% n = # rows/colums in distance matrix
% m = # parameters of the factor
% H = hat matrix
% I = I matrix
% G = Gower's centered matrix
% iter   = # iterations for permutation test

% -----Author:-----
% by David L. Jones, Aug-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% degrees of freedom:
df.among = m-1;
df.resid = n-m;
df.total = n-1;

% Sum-of-Squares:
SS.among = trace(H*G*H);
SS.resid = trace((I-H)*G*(I-H));
SS.total = trace(G);

% Mean Square:
MS.among = SS.among/df.among;
MS.resid = SS.resid/df.resid;

% pseudo-F:
F = MS.among/MS.resid;


%-----Permutation Tests:-----
if iter>0
   rand('state',sum(100*clock)); % set random generator to new state

   [nr,nc] = size(G);
   F_perm = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      G_perm = f_shuffle(G,2); % permute square symmetric matrix
            
      MS_perm  = trace(H*G_perm*H)/df.among;
      MSE_perm = trace((I-H)*G_perm*(I-H))/df.resid;
      F_perm(i) = MS_perm/MSE_perm;
   end;
   j = find(F_perm >= F);     % get randomized stats >= to observed statistic
   p = (length(j)+1)./(iter); % count values & convert to probability
else
   p = NaN;
end;
%-----------------------------

% wrap results up in a structure:
if (nargout > 5)
   result(1).so = {'Regression'};
   result(2).so = {'Residual'};
   result(3).so = {'Total'};
   
   result(1).df = df.among;
   result(2).df = df.resid;
   result(3).df = df.total;
   
   result(1).SS = SS.among;
   result(2).SS = SS.resid;
   result(3).SS = SS.total;
   
   result(1).MS = MS.among;
   result(2).MS = MS.resid;
   result(3).MS = NaN;
   
   result(1).F = F;
   result(2).F = NaN;
   result(3).F = NaN;
   
   result(1).p = p;
   result(2).p = NaN;
   result(3).p = NaN;
   
   result(1).R2 = SS.among/SS.total;
   result(2).R2 = SS.resid/SS.total;
   result(3).R2 = SS.total/SS.total;
   
end
