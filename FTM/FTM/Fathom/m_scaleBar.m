function m_scaleBar(lon,lat,km,wd,fnt,h,v)
% - create horizontal scale bar for M_Map plots
% 
% USAGE: m_scaleBar(lon,lat,km,wd,fnt,h,v)
% 
% lon = x coordinate of starting point of scale bar
% lat = y coodinate of starting point for scale bar
% km  = distance (km) between increments of scale bar
% wd  = width                                                (default = 3)
% fnt = font size                                            (default = 6)
% h   = proportion to narrow white bar for horizontal border (default = 0.45)
% v   = proportion to shorten white bar for vertical border  (default = 0.02)

% -----Notes:-----
% This function requires both the M_Map and Mapping Toolboxes and is typically
% called after producing an M_Map map.
% 
% This function can be tested by calculating the Earth's equatorial radius =
% 6,378.1370 km.
% 
% This function was written for North America, where longitudes are negative.

% -----References:-----
% http://en.wikipedia.org/wiki/Earth_radius

% -----Author:-----
% by David L. Jones, Apr-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2009: added support for wd, fnt, h, v

% -----Check input & set defaults:-----
if (nargin < 4), wd   = 3; end;  % set default width to 3
if (nargin < 5), fnt  = 6; end;  % set default font size to 6
if (nargin < 6), h  = 0.45; end; % default proportion to narrow width
if (nargin < 7), v  = 0.02; end; % set default proportion to shorten length
% -------------------------------------

% Get radius of the earth at the specified latitude:
a = 6378.1370; % equatorial radius (km)
b = 6356.7523; % polar radius (km)
r = sqrt( [(a^2 * cos(lat))^2 + (b^2 * sin(lat))^2] / ...
   [(a * cos(lat))^2 + (b * sin(lat))^2] );

% Specify how far apart to make segments:
deg = km2deg(km,r);

% Plot scale bar:
m_line([lon lon+(deg*2)],[lat lat],'color','k','LineWidth',wd);                   % black
m_line([lon+deg lon+deg+(deg*(1-v))],[lat lat],'color','w','LineWidth',wd*(1-h)); % white

% Add text labels:
m_text(lon, lat,sprintf('0\n\n'),'fontSize',fnt,'HorizontalAlignment', 'center', 'VerticalAlignment','middle');
m_text(lon+deg, lat, sprintf('%s\n\n', num2str(km)), 'fontSize',fnt,'HorizontalAlignment', 'center', 'VerticalAlignment','middle');
m_text(lon+(deg*2), lat, sprintf('   %s km\n\n', num2str(km*2)), 'fontSize',fnt,'HorizontalAlignment', 'center', 'VerticalAlignment','middle');
