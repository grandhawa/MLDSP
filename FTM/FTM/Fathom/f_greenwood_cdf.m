function Y = f_greenwood_cdf(G,n)
% - cumulative probability density function for Greenwood's statistic
%
% USAGE: Y = f_greenwood_cdf(G,n);
%
% G = Greenwood's statistic
% n = sample size
%
% Y = cumulative probability density
%
% SEE ALSO: f_greenwood, f_greenwood_pdf

% -----Notes:-----
% Algebraic simplification was performed using Wolfram Mathematica 9

% -----References:-----
% Burrows, P. M. 1979. Selected percentage points of Greenwood's statistic. J.
%  R. Statist. Soc. A 142(2): 256-258.
% Gardner, A. 1952. Greenwood's "Problem of Intervals": an exact solution for n
%  = 3. J. R. Statist. Soc. B 14: 135-139.
% Greenwood, M. 1946. The statistical study of infectious diseases. J. R.
%  Statist. Soc. A 109(2): 85-110.

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

% Replace out-of-bounds values with NaN's:
minG = 1/(n+1); % define minimum value of G
G(G<minG | G>1) = 0;

switch n
   case 1 % Burrows (1979):
      Y = sqrt(2*G-1);
      
   case 2 % Johnson (1946, attributed to L. Isserlis):
      % Y = 2/sqrt(3)*(G-(1/3)) .* (pi - 3*acos(1./(sqrt(6*G-2))) +...
      %   (3*sqrt(3)*sqrt(2*G-1))./(6*G-2));
      
      % Algebraic simplification:
      Y      = zeros(size(G));
      idx    = ~(G==minG);
      Y(idx) = (2*3^(1/2)*(G(idx) - 1/3).*(pi-3*acos(2^(1/2)./...
         (2*(3*G(idx)-1).^(1/2))) + (3*(6*G(idx) - 3).^(1/2))./(6*G(idx)-2)))/3;
      Y      = real(Y); % avoid complex numbers
      
   case 3 % Gardner (1952):
      Y = nan(size(G)); % initialize
      
      idx_1 = (G<=(1/3));
      if any(idx_1==1)
         % eq. P1:
         Y(idx_1) = 4*pi*(G(idx_1)-0.25).^(3/2);
      end
      
      idx_2 = (G>(1/3) & G<=(0.5));
      if any(idx_2==1)
         % Equation 7:
         Y(idx_2) = (2*pi)/sqrt(3) * ((3*G(idx_2)-(5/6)) -...
            2*(G(idx_2)-0.25).^(3/2)*sqrt(3));
      end
      
      idx_3 = (G>0.5);
      if any(idx_3==1)
         % Equation 5:
         % Y(idx_3) = 4*pi*(G(idx_3)-0.25).^(3/2)-8*pi*(G(idx_3)-0.25).^(3/2)...
         %   .*(1-((3*G(idx_3)-(5/6))./(4*(G(idx_3)-0.25).^(3/2)*sqrt(3))))...
         %   +12*(G(idx_3)-0.25).^(3/2).*(pi-2*asin((((2*G(idx_3)-0.5)./...
         %   (3*G(idx_3)-1))).^0.5)-(((3*G(idx_3)-(5/6)).*...
         %   acos(sqrt(6*G(idx_3)-2).^(-1)))./(2*(G(idx_3)-0.25).^(3/2)...
         %   *sqrt(3)))+((2*G(idx_3)-1).^(0.5)./(12*(G(idx_3)-0.25).^(3/2))));
         
         % Algebraic simplification:
         Y(idx_3) = (1/9)*(-27.206990463513264 - 56.548667764616276 *...
            sqrt(-0.25 + G(idx_3)) + 97.94516566864776 * G(idx_3) +...
            226.1946710584651 * sqrt(-0.25 + G(idx_3)) .* G(idx_3) +...
            9*(-1 + 2*G(idx_3)).^0.5 + 3*sqrt(3)*(5 - 18*G(idx_3)) .*...
            asec(sqrt(-2 + 6*G(idx_3))) - 216 * sqrt(-0.25 + G(idx_3)) .*...
            (-0.25 + 1*G(idx_3)) .* asin(((-0.5 + 2*G(idx_3)) ./...
            (-1 + 3*G(idx_3))).^0.5));
      end
   otherwise % Estimate probabilities using Johnson Curves:
      par = f_greenwood_par(n);                       % get parameters
      jsn = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
      Y   = f_johnson_cdf(G,jsn.coef,jsn.type);       % get probability
end

% Replace NaN's with 0's:
Y(isnan(Y)) = 0;
