function theta = f_vecAngle(a,b)
% - counter-clockwise angle between 2 points
%
% USAGE: theta = f_vecAngle(a,b);
%
% a = 2-d row vector for point a  (=  [Xa Ya])
% b = 2-d row vector for point b  (=  [Xb Yb])
%
% theta = angle between a & b in degrees
%
% See also: f_vecMagDir, f_vecTrans, f_vecUV, f_vecRotAngle

% ----- Notes: -----
% This program is used to determine the counter-clockwise angle
% (in degrees) between points A and B. The points are considered
% to be Cartesian coordinates of the heads of vectors starting at
% the origin.
%
% Note that the counter-clockwise angle from (A -> B) is not
% necessarily equal to that from (B -> A).
%
% This function is vectorized, so A and B may each be 2-d matrices
% specifying multiple pairs of points.

% ----- References: -----
% Feldman, M. 1997. The Win95 Game Programmer's Encyclopedia.
% Available from:
% http://www.geocities.com/SiliconValley/2151/win95gpe.html

% ----- Author(s): -----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: updated documentation

% ----- Check input: -----
if (size(a,2)~=2)
	error('A & B must each have 2 columns!');
end

if (size(a) ~= size(b))
	error ('A & B must be same size!');	
end
% ------------------------

% Convert to unit vectors:
[xA yA] = f_vecTrans(a(:,1),a(:,2),0,0,1,1);
[xB yB] = f_vecTrans(b(:,1),b(:,2),0,0,1,1);

% Find dir of each unit vector:
[magA dirA] = f_vecMagDir(xA,yA);
[magB dirB] = f_vecMagDir(xB,yB);

% Angle of rotation:
theta = dirB - dirA;
theta(find(theta<0)) = theta(find(theta<0)) + 360;

