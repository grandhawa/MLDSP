function mad = f_mad(X,c)
% - median absolute deviation
% 
% USAGE: mad = f_mad(X,c);
% 
% X = input matrix (rows = observations, cols = variables)
% c = adjustment factor (default = 1.4826)
% 
% mad = median absolute deviation
% 
% SEE ALSO: mad, f_outlier

% -----Notes:----
% This function is modelled after the R function 'mad', which uses a default
% adjustment factor to achiveve asymptotically normal consistency.

% -----Author:-----
% by David L. Jones, Aug-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), c =  1.4826; end % default adjustment factor

mad = c * median(abs(X-median(X)));
