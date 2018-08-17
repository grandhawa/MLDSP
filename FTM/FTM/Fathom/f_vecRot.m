function [uRot,vRot] = f_vecRot(u,v,theta,ccw,geo)
% - rotate vectors (U,V) by angle THETA
%
% USAGE: [uRot,vRot] = f_vecRot(u,v,theta,ccw,geo);
%
% u,v   = corresponding vector components
% theta = angle to vectors
% ccw   = rotate vectors counter-clockwise   (default = 1)
% geo   = use geographical coordinate system (default = 0)
%
% uRot,vRot = rotated vector components
%
% SEE ALSO: f_vecRotAngle, f_vecTrans, f_vecUV, f_vecMagDir

% ----- Notes: -----
% U,V components of velocity vectors can be extracted from data specifying
% only Magnitude and Direction using f_vecUV.
%
% THETA is an angle (in degrees) applied to rotate the vectors. When the
% vectors are rotated counter-clockwise (ccw) it is equivalent to rotating
% the coordinate system clockwise (and vise versa).
%
% It is useful to rotate the coordinate system clockwise, especially when 
% using a geographic coordinate system in order to align it with the local
% isobath coordinate system (i.e., the shoreline)---doing so partitions the
% vectors into cross-shore (u) and along-shore (v) components.
%
% If the geographic (= meteorological) coordinate system is used (GEO = 1),
% North = 0 degrees, East = 90 degress, and compass angles are measured in
% degrees clockwise from North. THETA can be obtained by taking the
% geographic compass angle (0=N,90=E) of the alongshore direction
% determined by the local isobaths using f_vecRotAngle.

% -----References:-----
% This code is based on F_VECPLOT, by the same author
%
% Formulas dealing with geographic coordinate systems from "Wind Direction
% Quick Reference", by Gordon Maclean; available from:
% http://www.atd.ucar.edu/rtf/facilities/isff/wind_ref.shtml

% ----- Author: -----
% by David L. Jones, Oct-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2005: added support for geographic coordinate system, clockwise and
%           counter-clockwise rotation.
% Jan-2005: updated documentation, changed '|' to '||'


% -----Check input & set defaults:-----
if (nargin < 4), ccw =  0; end; % default to ccw vector rotation
if (nargin < 5), geo =  0; end; % don't use geographical coordinates

if (size(u,1) ~= size(v,1))
   error('U and V must be same size!')
end

if (theta<0) || (theta>360)
   error('THETA must be positive angle b/n 0-360!');
end
% ---------------------------------------

nr = size(u,1); % # rows


if (geo==0)
   if (ccw==1)
      % Rotate from Geographic to Along-shore/Cross-shore coordinate system:
      % - vectors are rotated counter-clockwise
      % - coordinate system is rotated clockwise
      uRot = u * cos(theta*(pi/180)) - v * sin(theta*(pi/180));
      vRot = u * sin(theta*(pi/180)) + v * cos(theta*(pi/180));
   else
      % Rotate from Along-shore/Cross-shore to Geographic coordinate system:
      % - vectors are rotated clockwise
      % - coordinate system is rotated counter-clockwise
      uRot =  u * cos(theta*(pi/180)) + v * sin(theta*(pi/180));
      vRot = -u * sin(theta*(pi/180)) + v * cos(theta*(pi/180));
   end
else
   % Convert to mag/dir:
   [mag,dir] = f_vecMagDir(u,v);

   if (ccw==1)
      % Vectors are rotated counter-clockwise
      % Coordinate system is rotated clockwise
      dir = dir + theta;
      dir(find(dir>360)) = dir(find(dir>360)) - 360;
   else
      % Vectors are rotated clockwise
      % Coordinate system is rotated counter-clockwise
      dir = dir - theta;
      dir(find(dir<0))   = dir(find(dir<0)) + 360;
   end
      % Convert back to u,v:
      [uRot,vRot] = f_vecUV(mag,dir);
end


