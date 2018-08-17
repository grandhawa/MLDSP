function p = f_empPDF(X,Y,plt)
% - empirical probability density function
% 
% USAGE: p = f_empPDF(X,Y,plt);
% 
% X   = col vector of data forming empirical PDF
% Y   = col vector of data to evaluate
% plt = optionally plot empirical PDF  (default = 0)
% 
% p   = returns the probability Y is a member of the population represented
%       by the PDF
% 
% SEE ALSO: f_kde

% -----Author:-----
% by David L. Jones, Jan-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults:-----
if (nargin < 3), plt = 0;    end % default don't create plot

X  = X(:); % force column vector
Y  = Y(:); % force column vector

% Set interval for density estimates:
minVar   = min([X;Y]);
maxVar   = max([X;Y]);
rangeVar = maxVar - minVar;
minVar   = minVar - rangeVar/10;
maxVar   = maxVar + rangeVar/10;

% Use Botev's Kernel Density Estimator to generate probability densities:
[null,Di,Xi] = f_kde(X,[],minVar,maxVar,plt);

% Linearly interpolate position of Y along X's probability distribution:
p = interp1(Xi,Di,Y,'linear','extrap');
