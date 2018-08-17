function Y = f_greenwood_pdf(G,n)
% - probability density function for Greenwood's statistic
%
% USAGE: Y = f_greenwood_pdf(G,n);
%
% G = Greenwood's statistic
% n = sample size
%
% Y = probability density
%
% SEE ALSO: f_greenwood, f_greenwood_cdf

% -----Notes:-----
% Derivatives were calculated using Wolfram Mathematica 9

% -----References:-----
% Burrows, P. M. Selected percentage points of Greenwood's statistic. J. R.
%  Statist. Soc. A 142(2): 256-258
% Gardner, A. 1952. Greenwood's "Problem of Intervals": an exact solution for n
%  = 3. J. R. Statist. Soc. B 14: 135-139.

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
% Check size of input:
if (size(G,2)>1)
   error('G must be a column vector!');
end

% Check for missing values:
if any(isnan(G))
   error('G contains NaN''s!');
end

% Check N:
if ~isscalar(n), error('N must be a scalar!'); end
if (n-floor(n)>0), error('N must be a whole number'); end
if (n<0), error('N must be positive'); end
% -------------------------------------

minG = 1/(n+1); % get minimum value of statistic

% Replace out-of-bounds values with NaN's:
G(G< minG | G>1) = NaN;

switch n
   case 1
      Y      = zeros(size(G));         % initialize
      idx    = (G~=minG);              % prevent divide by 0
      Y(idx) = 1./sqrt(-1 + 2*G(idx)); % derivative of CDF
      
   case 2
      Y      = zeros(size(G));         % initialize
      idx    = (G~=minG);              % prevent divide by 0
      Y(idx) = (2*(pi-3*asec(sqrt(-2 + 6*G(idx)))))/sqrt(3); % derivative of CDF
      Y      = real(Y);                % avoid complex numbers
      Y(Y<0) = 0;                      % avoid negative probabilities
      
   case 3 % Gardner (1952)
      Y = nan(size(G)); % initialize
      
      idx_1 = (G<=(1/3));
      if any(idx_1==1)
         Y(idx_1) = 6*pi * sqrt(G(idx_1)-0.25);
      end
      
      idx_2 = (G>(1/3) & G<=(0.5));
      if any(idx_2==1)
         Y(idx_2) = 2*pi * (sqrt(3) - 3*(G(idx_2)-0.25).^0.5);
      end
      
      idx_3 = (G>0.5);
      if any(idx_3==1)
         Y(idx_3) = 2*pi*(sqrt(3) + 6*(G(idx_3)-0.25).^0.5) -...
            36*(G(idx_3)-0.25).^0.5 .* asin(((2*G(idx_3)-0.5)./...
            (3*G(idx_3)-1)).^0.5) - 6*sqrt(3)*(acos(((6*G(idx_3)-2)).^(-0.5))...
            + ((3*G(idx_3)-(5/6))./((6*G(idx_3)-3).^0.5.*(6*G(idx_3)-2)))) +...
            ((18*G(idx_3)-5) ./ ((2*G(idx_3)-1).^0.5.*(6*G(idx_3)-2)));
      end
   otherwise % Estimate probabilities using Johnson Curves:
      par = f_greenwood_par(n);                       % get parameters
      jsn = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
      Y   = f_johnson_pdf(G,jsn.coef,jsn.type);       % get probability
end
