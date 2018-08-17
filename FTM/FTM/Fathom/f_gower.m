function G = f_gower(dis)
% - Gower's centered matrix
% 
% USAGE: G = f_gower(dis)
% 
% dis = symmetric distance matrix
% G   = Gower's centered matrix
% 
% SEE ALSO: f_dis

% -----References:-----
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.

% -----Author:-----
% by David L. Jones, Oct-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (f_issymdis(dis) == 0)
   error('Input DIS must be a square symmetric distance matrix!');
end
% ----------------------

n = size(dis,1); % # rows

% Centering terms:
I   = eye(n,n);
uno = ones(n,1);

A = -0.5*(dis.^2);
G = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno')); % Gower's centered matrix
