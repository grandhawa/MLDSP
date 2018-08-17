function [p,p_labels] = f_poly(x,y,d)
% - create a matrix of polynomials of degree d
%
% USAGE: [p,p_labels] = f_poly(x,y,degree)
%
% x,y = column vectors
% d   = degree of polynomial (= 2 or 3)
%
% p        = matrix of polynomials
% p_labels = cell array of labels

% -----Author:-----
% by David L. Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

x = x(:); % force column vectors
y = y(:);

n = size(x,1);

if (n ~= size(y,1))
   error('X & Y must have same # of rows!');
end

if (d~=3) & (d~=2)
   error('Only polynomials of degree 2 or 3 are supported');
end

if (d==3) % cubic polynomial
   p = zeros(n,9);
   
   p(:,1) = x;
   p(:,2) = y;
   p(:,3) = x.^2;
   p(:,4) = x.*y;
   p(:,5) = y.^2;
   p(:,6) = x.^3;
   p(:,7) = (x.^2).*y;
   p(:,8) = x.*(y.^2);
   p(:,9) = y.^3;
   
   p_labels = {'X' 'Y' 'X^2' 'XY' 'Y^2' 'X^3' 'X^2Y' 'XY^2' 'Y^3'};
   
else % (d==2) quadratic polynomial
   p = zeros(n,5);
   
   p(:,1) = x;
   p(:,2) = y;
   p(:,3) = x.^2;
   p(:,4) = x.*y;
   p(:,5) = y.^2; 
   
   p_labels = {'X' 'Y' 'X^2' 'XY' 'Y^2'};
end
   
