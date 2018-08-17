function [tx,ty] = f_vecTrans(x,y,rot,transl,scale,unit)
% - transform 2d vector coordinates
%
% USAGE: [tx,ty] = f_vecTrans(x,y,rot,transl,scale,unit);
%
% x,y    = column vectors specifying coordinate pairs
% rot    = angle of rotation in degrees   (default = 0)
% transl = translation    [dx dy]         (default = [0 0])
% scale  = scaling factor [sx sy]         (default = [1 1])
% unit   = convert vectors to unit length (default = 0)
%
% tx,ty  = transformed coordinates
%
% See also: f_vecAngle, f_vecMagDir, f_vecTrans3d, f_vecUV

% ----- Notes:-----
% This function is used to rotate, translate, and/or scale 2
% dimensional Cartesian coordinates. If the input coordinates
% specify the heads of (velocity) vectors, they may additionally
% be converted to unit length.
%
% X and Y coordinates can be translated or scaled asymmetrically
% if 2 values are specified for these parameters. If only 1 value
% is provided, coordinates are translated or scaled symmetrically.

% ----- References: -----
% Feldman, M. 1997. The Win95 Game Programmer's Encyclopedia.
% Available from:
% http://www.geocities.com/SiliconValley/2151/win95gpe.html

% ----- Author(s): -----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% ----- Check input & set defaults: -----
if (nargin < 3), rot    = 0;     end; % no rotation
if (nargin < 4), transl = [0 0]; end; % no translation
if (nargin < 5), scale  = [1 1]; end; % no scaling
if (nargin < 6), unit   = 0;     end; % no unit length transformation

if (size(x,2)~=1)
	error('X & Y must be column vectors!');
end

if (size(x) ~= size(y))
	error ('X & Y must be same size!');	
end

if (rot<0) | (rot>360)
	error('Angle of rotation must be between 0-360!')
end

if (sum(find(scale==0))>0)
	error('You can''t have a scaling factor of 0!');
end

transl = transl(:)'; % make sure it's a row vector
if (size(transl,2)==1)
	transl = [transl transl]; % symmetric translation
end

scale = scale(:)'; % make sure it's a row vector
if (size(scale,2)==1)
	scale = [scale scale];    % symmetric scaling
end
% ---------------------------------------

[nr,nc] = size(x);

% Convert to unit vector:
if (unit>0)
	vecLength = sqrt(x.^2 + y.^2); % get vector length
	x = x .* (vecLength.^-1);      % divide each component by length
	y = y .* (vecLength.^-1); 
end

% Initialize as identity matrices (no effect):
I = eye(3,3); % identity matrix
rMat = I;
tMat = I;
sMat = I;

% Rotation matrix:
if (rot ~= 0)
	rot = rot * pi/180; % degrees to radians
	rMat(1,1) =  cos(rot);
	rMat(1,2) = -sin(rot);
	rMat(2,1) =  sin(rot);
	rMat(2,2) =  cos(rot);
end

% Translation matrix:
tMat(1,3) = transl(1);
tMat(2,3) = transl(2);

% Scaling matrix:
sMat(1,1) = scale(1);
sMat(2,2) = scale(2);

% Transformation matrix:
transMat = I * rMat * tMat * sMat;

% Apply transformation:
xy = [x y ones(nr,1)]; % add column of 1's
txy = (transMat * xy')';

% Extract output:
tx = txy(:,1);
ty = txy(:,2);

