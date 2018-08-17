function [mag,dir] = f_vecMagDir(u,v,geo)
% - get magnitude & direction from u,v vector components
%
% Usage: [mag,dir] = f_vecMagDir(u,v,geo);
%
% u,v = column vectors of Cartesian coordinates of heads
%       of vector components
% geo = use geographical coordinate system (default = 0)
%
% mag = length of vector
% dir = angle of rotation (in degrees between 0-360)
%
% See also: f_vecAngle, f_vecTrans, f_vecUV

% ----- Notes: -----
% This function is used to obtain Polar coordinates (magnitude
% & direction) of a vector given its Cartesian coordinates (U & V
% vector components). The direction is the counter-clockwise angle
% of rotation.
%
% The programs uses the Matlab function ATAN2 which relies on the
% sign of both input arguments to determine the quadrant of the result.
%
% If the geographic (= meteorological) coordinate system is used, North =
% 0 degrees, East = 90 degress, and compass angles are measured in degrees
% clockwise from North.

% -----Reference:-----
% Formulas dealing with geographic coordinate systems from "Wind Direction
% Quick Reference", by Gordon Maclean; available from:
% http://www.atd.ucar.edu/rtf/facilities/isff/wind_ref.shtml

% ----- Author(s): -----
% by David L. Jones, Sept-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2002: tune up, check input, added notes, vectorized code
% Oct-2005: included support for geographic coordinate system

% -----Check input & set defaults:-----
if (nargin < 3), geo =  0; end; % don't use geographical coordinates

u = u(:); % make sure they're column vectors
v = v(:);

if (size(u,1) ~= size(v,1))
	error('U and V must be of same size!');
end
% ----------------------

if (geo == 0)
   mag = sqrt(u.^2 + v.^2);
   dir = (atan2(v,u)) .* 180/pi;
   dir(find(dir>360)) = dir(find(dir>360)) - 360;
   dir(find(dir<0))   = dir(find(dir<0)) + 360;
else
   mag = sqrt(u.^2 + v.^2);
   dir = (atan2(u,v)) .* 180/pi;
   dir(find(dir>360)) = dir(find(dir>360)) - 360;
   dir(find(dir<0))   = dir(find(dir<0)) + 360;
end

% return column vectors
mag = mag(:);
dir = dir(:);
