function [e,n] = f_ll2utm(lon,lat)
% - legacy function (use f_deg2utm instead)
%
% USAGE: [e,n] = f_ll2utm(lon,lat);
%
% lon = longitude in decimal degrees (column vector)
% lat = latitude in decimal degrees  (column vector)
%
% e = easting
% n = northing
%
% This function requires the M_Map Toolbox
% 
% SEE ALSO: f_deg2utm, f_utm2deg

% -----Notes:-----
% This is a quick and dirty approach utilizing the M_Map Toolbox to
% generate UTM coordinates. This was done primarily to obtain geographic
% coordinates that maintain a constant distance relationship. This is
% necessary if one is interested in using the coordinates in a spatial
% analysis (e.g., lon/lat as (co)-variables in RDA, multiple linear
% regression, etc.).
%
% This function will return UTM-like coordinates, but their actual values
% have not been centered according to any actual UTM zone. This makes
% little difference, though, as these variables are usually centered or
% standardized during RDA, etc.

% -----Author:-----
% by David L. Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
lon = lon(:); % force column vector
lat = lat(:);

if (size(lon,1) ~= size(lat,1))
   error('LON and LAT must have same # rows!')
end

% Create a base map:
m_proj('UTM','longitudes',[min(lon) max(lon)],'latitudes',[min(lat) max(lat)]);

m_grid('box','normal','fontsize',8,'linestyle','none');

% Plot points:
sites_hdl = m_line(lon,lat,'marker','o','markersize',4,...
   'MarkerFaceColor','r','MarkerEdgeColor','none','LineStyle','none');

daspect([1 1 1]);

% Extract UTM-like coordinates:
e = get(sites_hdl,'XData');
n = get(sites_hdl,'YData');

e = e(:); % force column vector
n = n(:);

