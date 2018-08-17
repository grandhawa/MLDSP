function par = f_greenwood_par(n)
% - get parameters of the distribution of Greenwood's statistic for n>2
% 
% USAGE: result = f_greenwood_par(n)
% 
% n   = sample size
% par = parameters of the distribution as: [mean sd skewness kurtosis]

% -----Notes:-----
% SKEW: 3rd central moment divided by the cube of the SD.
% KURT: 4th central moment divided by fourth power of its SD.

% -----References:-----
% Moran, M. J. 1947. The random division of an interval. J. R. Statisti. Soc. B
% 9: 92?98. [Corrigendum. 1981. J. R. Statisti. Soc. A 144: 388].

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if any(n<3), error('N must be 3 or greater!'); end
if any(n-floor(n)>0), error('N must be a whole number'); end
% -------------------------------------

n = n(:); % force as row vector

% Get the 1st moment and 2nd-4th 'central' moments: (Moran, 1947 p.94):
u1 = 2*(n+2).^(-1);
u2 = 2^2*n.*(n+2).^(-2).*((n+3).*(n+4)).^(-1);
u3 = 2^3*(10*n.^2-4*n).*(n+2).^(-3).*((n+3).*(n+4).*(n+5).*(n+6)).^(-1);
u4 = 2^4*(3*n.^4 + 303*n.^3 + 42*n.^2 - 24*n).*(n+2).^(-4).*((n+3).*(n+4).*...
   (n+5).*(n+6).*(n+7).*(n+8)).^(-1);

% Get parameters of the distribution:
mu = u1;          % mean
sd = sqrt(u2);    % standard deviation
sk = u3./(sd).^3; % kurtosis
ku = u4./(sd).^4; % skewness

% Format results for output:
par = [mu sd sk ku];
