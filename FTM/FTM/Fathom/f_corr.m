function [r,p] = f_corr(x,y,rank,iter,tail)
% - Pearson's, Spearman's, or Kendall's correlation between 2 vectors
%
% USAGE: [r,p] = f_corr(x,y,rank,iter,tail);
%
% x,y  = input vectors
% rank = rank correlation
%           (default = 0, Pearson)
%           (1 = Spearman, 2 = Kendall)
%
% iter = iterations for permutation test (default = 0)
% tail = type of permutation test
%           (1 = one-tailed)
%           (2 = two-tailed, default)
%
% r = correlation coefficient
% p = randomized significance of r
% 
% SEE ALSO: corrcoef

% This function is used to calculate 'r', the correlation coefficient
% between 2 input vectors with an optional permutation-based
% significance test. The significance of 'r' is evaluated through
% permutation of the elements in X.

% -----Author:-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam. (pp.553)
%
% Sokal, R. R. and F. J. Rohlf. 1995. Biometry - The principles and practice
% of statistics in bioligical research. 3rd ed. W. H. Freeman, New York.
% xix + 887 pp. (pp. 594-595 for Kendall's tau)

% 02-Apr-2002: added Kendall's tau; 1- & 2-tailed tests
% Apr-2003: edited documentation
% Mar-2010: edited documentation
% Mar-2011: f_transform no longer needs to be transposed

% -----Set defaults & check input:-----
if (nargin < 3), rank = 0; end;          % don't rank by default
if (nargin < 4), iter = 0; p = NaN; end; % don't perform permutation test
if (nargin < 5), tail = 2; end;          % 2-tailed test by default

if (sum([0 1 2] == rank)<1); error('Invalid option for RANK'); end;
if (sum([1 2]   == tail)<1); error('Invalid option for TAIL'); end;

% make sure they're column vectors
x = x(:);
y = y(:);

n = length(x);

if (n ~= length(y))
   error('Input vectors must same size!');
end
% -------------------------------------

% Optionally rank:
if (rank>0)
   x = f_ranks(x);
   y = f_ranks(y);
end

if (rank<2) %%--Pearson's & Spearman's R--%%

   x = f_transform(x,7); % standardize
   y = f_transform(y,7);

   % Take sum of cross-products, divide by n-1:
   r = (sum(x.*y))/(n-1);

else %%--Kendall's tau--%%

   xy = sortrows([x y],1);    % sort according to x
   x  = xy(:,1); y = xy(:,2); % extract sorted x & y

   c = 0; % initialize counter
   for i = 1:(n-1)
      % count all subsequent ranks > than the one being considered
      c = c + sum(y(i+1:end)>y(i));
      c = c + sum(0.5*(y(i+1:end)==y(i))); % ties count as only half
   end;

   NN = 4*c - (n*(n-1));

   % # ties in x & y:
   xTies = sub_ties(x);
   yTies = sub_ties(y);

   r = NN/(sqrt((n*(n-1) - xTies) * (n*(n-1) - yTies))); % tau
end

%-----Permutation Test:-----
if iter>0
   randStat = zeros(iter-1,1); % preallocate results array
   for i = 1:(iter-1) % iter-1 since observed statistic is an iteration
      xx = (f_shuffle(x)); % permute x
      if (rank<2) % Pearson's or Spearman's R
         randStat(i) = f_corr(xx,y,rank,0); % collect randomized stat
      else % Kendall's tau
         randStat(i) = f_corr(xx,y,2,0);    % collect randomized stat
      end
   end
   if (tail>1)
      j = find(randStat.^2 >= r^2); % 2-tailed test of sign. correlation
   elseif (r<0)
      j = find(randStat <= r); % 1-tailed test of sign - correlation
   else
      j = find(randStat >= r); % 1-tailed test of sign + correlation
   end
   p = (length(j)+1)./(iter);  % count vales & convert to probability
else
   p = NaN;
end
%-----------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        SUBFUNCTION           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noTied = sub_ties(z)
% - returns the number of ties in z
zz  = unique(z);
noZ = length(zz);
noTied = 0; % initialize variable
for i = 1:noZ
   noOccurred = length(find(z == zz(i)));
   if noOccurred>1,
      noTied = noTied + noOccurred;
   end
end



