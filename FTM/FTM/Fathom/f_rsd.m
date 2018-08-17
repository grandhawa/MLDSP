function rsd = f_rsd(X)
% - percent relative standard deviation
% 
% USAGE: rsd = f_rsd(X)
% 
% x   = input matrix (rows = observations, columns = variables)
% rsd = row vector with values for 'percent relative standard deviation' for
%       each column (variable) in X
% 
% SEE ALSO f_prd, f_perError

% -----NOTES:-----
% In analytical chemistry 'Percent RDS' is used to measure the precision
% (= repeatability) of an assay. 

% -----Author:-----
% by David L. Jones, Aug-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

rsd = (std(X)./mean(X))*100;
