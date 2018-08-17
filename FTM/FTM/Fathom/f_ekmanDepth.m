function z = f_ekmanDepth(lat,ws)
% - depth of the Ekman layer
%
% USAGE: z = f_ekmanDepth(lat,ws)
%
% lat = latitude
% ws  = wind speed in m/s
% z   = depth in meters
%
% See also: f_windStress

% ----- Notes: -----
% This function is used to calculate the depth of
% the Ekman layer as a function of wind speed and
% latitude.

% ----- Author(s): -----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

lat = lat * pi/180; % degrees to radians
z   = (4.3*(ws)) /sqrt(sin(lat));
