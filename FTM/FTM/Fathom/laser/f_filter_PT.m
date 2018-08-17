function Xf = f_filter_PT(X)
% - filter/smooth 'otolith profile' surface transect data
%
% USAGE: Xf = f_filter_PT(X)
%
% X  = structure of 'otolith profile' surface transect data originally created
%      by f_cps2ppm_PT
%
% Xf = filtered/smoothed version of X
%
% SEE ALSO: f_cps2ppm_PT, f_extract_PT, f_filterSinclair

% -----Notes:-----
% This function applies the 'f_filterSinclair' function to the data in the PPM,
% LOD, and RATIO fields in a structure originally created by the f_cps2ppm_PT
% function. It will also work with structures processed with f_extract_PT
% function, etc.

% -----Author:-----
% by David L. Jones, Sep-2013
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
% Check whether data have been filtered/smoothed:
if isequal(X.sm,{'data filtered'})
   error('Check SM: data have already been filtered!')
end
% ----------------------

% Filter/smooth times-series data:
X.ppm   = f_filterSinclair(X.ppm);
X.LOD   = f_filterSinclair(X.LOD);
X.ratio = f_filterSinclair(X.ratio);

% Rename structure for output:
f_rename X Xf;
