function [u,v] = f_vecUV(mag,dir,geo)
% - returns U,V components of a vector given its magnitude & direction
%
% Usage: [u,v] = f_vecUV(mag,dir,geo);
%
% mag = column vector specifying magnitude of vectors (in arbitrary units)
% dir = column vector indicating angle of rotation (in degrees from 0-360)
% geo = use geographic coordinate system (default = 0)
%
% u,v = Cartesian coordinates of heads of vector components
%
% See also: f_vecAngle, f_vecMagDir, f_vecTrans

% -----References:-----
% Formulas dealing with geographic coordinate systems from "Wind Direction
% Quick Reference", by Gordon Maclean; available from:
% http://www.atd.ucar.edu/rtf/facilities/isff/wind_ref.shtml

% -----Author:-----
% by David L. Jones, Sept-2001
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Notes:-----
% This function is used to obtain Cartesian coordinates (the U,V vector
% components) of a vector given its Polar coordinates (magnitude &
% direction).
%
% If the geographic (= meteorological) coordinate system is used, North =
% 0 degrees, East = 90 degress, and compass angles are measured in degrees
% clockwise from North.

% Nov-2002: tune-up, check input, added notes
% Oct-2005: included support for geographic coordinate system
% Oct-2009: changed '|' to '||'

% -----Check input & set defaults:-----
if (nargin < 3), geo =  0; end; % don't use geographical coordinates

mag = mag(:); % make sure they're column vectors
dir = dir(:);

if (size(mag,1) ~= size(dir,1))
	error('MAG and DIR must be of same size!');
end

if (min(dir)<0) || (max(dir>360))
	error('DIR must range from 0-360 degrees');
end
% ----------------------

if (geo == 0)
   u = mag .* cos(dir * pi/180);
   v = mag .* sin(dir * pi/180);
else
  u = mag .* sin(dir * pi/180);
  v = mag .* cos(dir * pi/180);
end

