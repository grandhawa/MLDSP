function cmap = f_colormap(n)
% - create a custom colormap
% 
% USAGE: cmap = f_colormap(n);
% 
% n    = # color levels
% cmap = matrix of colors to use with colormap command
 
% -----Notes:-----
% This preliminary version only does shades of gray, but varies the values
% for rgb to colors; 

% -----Author:-----
% by David L. Jones,<djones@rsmas.miami.edu> Aug-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

c    = [0:1/n:1]';
cmap = repmat(c,1,3);
