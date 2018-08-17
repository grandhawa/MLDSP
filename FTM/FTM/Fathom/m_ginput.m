function crds = m_ginput
% - get lon/lat coordinates from an M_Map figure
%
% USAGE: crds = m_ginput;
%
% crds = [lon lat] coordinates
%
% SEE ALSO: ginput

% -----Notes:-----
% This program is used to provide 'ginput' functionality
% when working with a figure created with the M_Map Toolbox.
% 
% It requires the M_Map Toolbox.

% -----Author:-----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

xy        = ginput;
[lon,lat] = m_xy2ll(xy(:,1),xy(:,2));
crds      = [lon lat];
