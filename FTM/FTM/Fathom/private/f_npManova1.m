function [df,SS,MS,F,p,result] = f_npManova1(n,m,H,I,G,iter,verb)
% - utility function called by f_npManova (1-way MANOVA)
% 
% n    = # rows/colums in distance matrix
% m    = # parameters of the factor
% H    = hat matrix
% I    = I matrix
% G    = Gower's centered matrix
% iter = # iterations for permutation test
% verb = verbose output

% -----Author:-----
% by David L. Jones
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2008: replaced sum(diag()) with trace
% Mar-2010: added verb option to minimize output from f_npManovaPW
% Dec-2010: replaced 'for' with 'parfor' in permutation test
% Nov-2011: added code for plotting distribution of permuted F.
% Mar-2013: replaced 'parfor' with 'for' in permutation test to avoid error
%           in new Java

% Degrees of freedom:
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

% Pseudo-F:
F = MS.among/MS.resid;

%-----Permutation tests:-----
if (iter>0)
   if (verb>0)
      fprintf('\nPermuting the data %d times...\n',iter-1);
   end
   % [nr,nc] = size(G);
   randStat = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      Gpermed     = f_shuffle(G,2); % permute square symmetric matrix
      among       = trace(H*Gpermed*H)/df.among;
      resid       = trace((I-H)*Gpermed*(I-H))/df.resid;
      randStat(i) = among/resid;
   end
   j = find(randStat >= F);   % get randomized stats >= to observed statistic
   p = (length(j)+1)./(iter); % count values & convert to probability
else
   p = NaN;
end
%-----------------------------

% % -----Plot the randomized distribution:-----
% size(randStat)
% size(F)
% 
% 
% dave = [randStat;F];
% figure;
% hist(dave,100)
% xlabel('F');
% ylabel('Frequency')
% % ---------------------------------
 
% Wrap results up in a structure for 1-way MANOVA's:
if (nargout > 5)
   result(1).so = {'factor 1'};
   result(2).so = {'residual'};
   result(3).so = {'total'};
   
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
end
