function angle = f_vecRotAngle(data)
% - determine the angle of an isobath
%
% USAGE angle = f_vecRotAngle(data)
%
% data  = coordinates returned by ginput after clicking on 2 points along an
%         isobath contour created by, say, M_Map; the first point should be
%         DOWNSTREAM of the 2nd point
%
% angle = degrees to rotate vectors counter-clockwise to align with local
%         isobath system
%
% SEE ALSO: f_vecRot, f_vecAngle

% ----- Author(s): -----
% by David L. Jones, Nov-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Dec-2009: edited for more generic use

% translate point 2 to it is a vector coming from the origin:
pnt2 = [data(2,:) - data(1,:)];

% make point 1 vector of, say, mag=10, dir=90:
pnt1 = [0 10];

% clockwise angle of rotation from 1 -> 2:
angle = f_vecAngle(pnt1,pnt2);
