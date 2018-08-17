function c = m_extContours(lon,lat,depth,interval)
% - extract contours from an M_Map map
%
% USAGE: c = m_extContours(lon,lat,depth,interval);
%
% lon,lat  = coordinates of extracted elevation/bathymetry data to contour
% depth    = corresponding depths
% interval = [min max] range of extracted contours (e.g., [-10 -10]
%
% c = [lon lat] coordinates of extracted contour
% 
% SEE ALSO: readmeUSCRM.m

% -----Notes:-----
% This program requires the M_Map Toolbox

% -----Author:-----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2009: almost a complete re-write, no longer querying 'Xdata' using the
%           'get' command

% Set up projection:
m_proj('mercator','longitudes',[min(lon) max(lon)],'latitudes',[min(lat) max(lat)]);

% Create contour plot:
[cs,h] = m_contour(lon,lat,depth,interval);

cs      = cs';
idx     = find( (abs(cs)>=1)==1 ); % eliminate contour headers
cs(idx) = NaN;

plot(cs(:,1),cs(:,2),'b-');

% Convert x,y to lon,lat:
[lon,lat] = m_xy2ll(cs(:,1),cs(:,2));

c = [lon lat];
