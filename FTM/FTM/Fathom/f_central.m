function [result,se] = f_central(x,method,zed,arc,err)
% - central tendency (mean, median, or mode)
%
% USAGE: [result,se] = f_central(x,method,zed,arc,err);
%
% x      = input data (rows = objects, cols = variables)
% method = mean (= 1), median(= 2), or mode (= 3)
% zed    = include zeros (default = 1)
% arc    = transform the result via the arcsine (= 1) or none (= 0, default)
% err    = calculate standard error (default = 0)
% 
% result = measure of central tendency
% se     = standard error

% -----Notes:-----
% Support for not including 0's when calculating the central tendency was added
% in order to use the so-called 'delta approach' for calculating means, etc.
% of 'abundance when present'
% 
% Support for arcsine transformation was added in order to return an arcsine
% transformed mean of presence/absence data when working with frequencies or
% proportions. For example: f_central(f_normal(x,'01'),1,1,1)

% -----Author:-----
% by David L. Jones, Oct-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2008: included support for deleting 0's in input
% Jan-2009: support for arcsine transformation
% Oct-2009: support for standard error

% -----Check input & set defaults:-----
if (nargin < 3), zed = 1; end % include 0's by default
if (nargin < 4), arc = 0; end % no arcsine transformation by default
if (nargin < 5), err = 0; end % no standard error by default

% Remove 0's before calculating statistics:
if (zed==0)
   x(x==0) = [];
end

% Calculate statistic:
switch method
   case 1
      result = mean(x);

   case 2
      result = median(x);
      
   case 3 
      result = mode(x);
      
   otherwise
      error('METHOD must be 1, 2, or 3')
end

% Arcsine transform:
if (arc>0)
   result = f_normal(result,'ang');
end


% Standard error:
if (err>0)
   se = f_stdErr(x);
else
   se = NaN;
end
