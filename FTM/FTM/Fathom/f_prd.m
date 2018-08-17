function prd = f_prd(X,Y)
% - percent relative difference between 2 sets of values
% 
% USAGE: prd = f_prd(X,Y)
% 
% X = vector or matrix of one set of values
% Y = vector or matrix of corresponding set of values to compare
% 
% prd = percent relative difference
% 
% SEE ALSO f_perError, f_rsd

% -----NOTES:-----
% 'Percent Difference' is used to measure the degree of divergence of 2 values
% expected to be numerically equivalent, such as repeated measurements of the
% same thing or calculations/estimates based on different methods, equations,
% etc.

% -----References:-----
% http://en.wikipedia.org/wiki/Percent_difference

% -----Author:-----
% by David L. Jones, Aug-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if size(X,1) ~= size(Y,1)
   error('X and Y must have same number of rows!');
end

if size(X,2) ~= size(Y,2)
   error('X and Y must have same number of columns!');
end
% ----------------------

prd = abs((X - Y)./((X+Y)/2)) *100;
