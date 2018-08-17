function seaIce = f_NSIDCinterp(result,lon,lat)
% - interpolate sea ice coverage values for a specific location
% 
% USAGE: seaIce = f_NSIDCinterp(result,lon,lat)
% 
% result = structure created by the f_NSIDCimport function
% lon    = column vector of longitudes to interpolate sea ice values 
% lat    = column vector of latitudes  to interpolate sea ice values 
% 
% seaIce = sea ice coverage (0-100%) corresponding to lon/lat coordinates
% 
% SEE ALSO: f_NSIDCimport

% -----Author:-----
% by David L. Jones, July-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Force column vectors:
lon = lon(:);
lat = lat(:);

% Unwrap input:
X = result.X(:);
Y = result.Y(:);
P = result.P(:);

% Set map projection:
m_proj('stereographic','long',0,'lat',-90,'radius',52);

% Convert lon,lat coordinates to map coordinates:
[mapX,mapY]=m_ll2xy(lon,lat,'clip','off');

% % Construct the interpolant:
F = TriScatteredInterp(X,Y,P);

% Interpolate Sea Ice Coverage values:
seaIce = F(mapX,mapY);
