function pe = f_perError(X,Y)
% - percent error between measured and accepted values
% 
% USAGE: pe = f_perError(X,Y)
% 
% X  = vector or matrix of accepted values 
% Y  = vector or matrix of corresponding measured values
% 
% pe = percent error
% 
% SEE ALSO f_prd, f_rsd

% -----NOTES:-----
% 'Percent Error' is used to measure the accuracy of an assay, by comparing an
% experimental (measured) value with a theoretical (accepted) value.

% -----References:-----
% http://en.wikipedia.org/wiki/Percent_difference

% -----Author:-----
% by David L. Jones, Sep-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if size(X,1) ~= size(Y,1)
   error('Y and X must have same number of rows!');
end

if size(X,2) ~= size(Y,2)
   error('Y and X must have same number of columns!');
end
% ----------------------

pe = ((Y - X)./X)*100;
