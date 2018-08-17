function result = f_cov(x,method)
% - returns covariance (= dispersion) matrix from X
%
% USAGE: result = f_cov(x,method)
%
% x      = input matrix
% method = divide by degrees of freedom (= 1, default) 
%
% SEE ALSO: cov

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam.

% -----Author(s):-----
% by David L. Jones, Sep-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

if (nargin < 2), method = 1; end; % divide by degrees of freedom


n = size(x,1);

% Center data on column mean:
Xctr = f_center(x);

if (method>0)
   result = 1/(n-1)*(Xctr'*Xctr); % covariance matrix (eq. 4.6 page 392)
else
   result = Xctr'*Xctr;
end