function [tx,ty,tz] = f_vecTrans3d(x,y,z,xrot,yrot,zrot,transl,scale,unit);
% - transform 3d vector coordinates
%
% USAGE: [tx,ty,tz] = f_vecTrans3d(x,y,z,xrot,yrot,zrot,transl,scale,unit);
%
% x,y,z  = column vectors specifying coordinate triplets
% xrot   = rotation about X-axis in degrees (default = 0)
% yrot   = rotation about Y-axis in degrees (default = 0)
% zrot   = rotation about Z-axis in degrees (default = 0)
% transl = translation    [dx dy dz]        (default = [0 0 0])
% scale  = scaling factor [sx sy sz]        (default = [1 1 1])
% unit   = convert vectors to unit length   (default = 0)
%
% tx,ty,tz  = transformed coordinates
%
% See also: f_vecAngle, f_vecMagDir, f_vecTrans, f_vecUV, f_vecPlot

% ----- Notes:-----
% This function is used to rotate, translate, and/or scale 3
% dimensional Cartesian coordinates. If the input coordinates
% specify the heads of (velocity) vectors, they may additionally
% be converted to unit length.
%
% X, Y, and Z coordinates can be translated or scaled asymmetrically
% if 3 values are specified for these parameters. If only 1 value
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
if (nargin < 3), xrot   = 0;       end; % no x rotation
if (nargin < 4), yrot   = 0;       end; % no y rotation
if (nargin < 5), zrot   = 0;       end; % no z rotation
if (nargin < 6), transl = [0 0 0]; end; % no translation
if (nargin < 7), scale  = [1 1 1]; end; % no scaling
if (nargin < 8), unit   = 0;       end; % no unit length transformation

if (size(x,2)~=1)
	error('X, Y, & Z must be column vectors!');
end

if (size(x) ~= size(y)) | (size(x) ~= size(z))
	error ('X, Y, & Z must be same size!');	
end

if (sum(find([xrot yrot zrot]<0))>0) | (sum(find([xrot yrot zrot]>360))>0)
	error('Angle of rotation must be between 0-360!')
end

if (sum(find(scale==0))>0)
	error('You can''t have a 0 scaling factor!');
end

transl = transl(:)'; % make sure it's a row vector
switch size(transl,2)
case 1
	transl = [transl transl transl]; % symmetric translation
case 2
	error('You only specified translation for X & Y');
case 3
	% all 3 specified
otherwise
	error('You specified >3 values for X,Y,Z translation')
end

scale = scale(:)'; % make sure it's a row vector
switch size(scale,2)
case 1
	scale = [scale scale scale]; % symmetric translation
case 2
	error('You only specified scaling for X & Y');
case 3
	% all 3 specified
otherwise
	error('You specified >3 values for X,Y,Z scaling!')
end
% ---------------------------------------

[nr,nc] = size(x);

% Convert to unit vector:
if (unit>0)
	vecLength = sqrt(x.^2 + y.^2 + z.^2); % get vector length
	% divide each component by length:
	x = x .* (vecLength.^-1);             
	y = y .* (vecLength.^-1); 
	z = z .* (vecLength.^-1); 
end

% Initialize as identity matrices (no effect):
I = eye(4,4); % identity matrix
xMat = I;
yMat = I;
zMat = I;
tMat = I;
sMat = I;

% ----- Rotation Matrices: -----

% About the X-Axis:
if (xrot ~= 0)
	xrot      = xrot * pi/180; % degrees to radians
	xMat(2,2) =  cos(xrot);
	xMat(2,3) = -sin(xrot);
	xMat(3,2) =  sin(xrot);
	xMat(3,3) =  cos(xrot);
end

% About the Y-Axis:
if (yrot ~= 0)
	yrot     = yrot * pi/180;
	yMat(1,1) =  cos(yrot);
	yMat(1,3) =  sin(yrot);
	yMat(3,1) = -sin(yrot);
	yMat(3,3) =  cos(yrot);
end

% % About the Z-Axis:
if (zrot ~= 0)
	zrot      = zrot * pi/180;
	zMat(1,1) =  cos(zrot);
	zMat(1,2) = -sin(zrot);
	zMat(2,1) =  sin(zrot);
	zMat(2,2) =  cos(zrot);
end

% -----Translation Matrix: -----
tMat(1,4) = transl(1);
tMat(2,4) = transl(2);
tMat(3,4) = transl(3);


% ----- Scaling Matrix: -----
sMat(1,1) = scale(1);
sMat(2,2) = scale(2);
sMat(3,3) = scale(3);

% Transformation matrix:
transMat = I * xMat * yMat * zMat* tMat * sMat;

% Apply transformation:
xyz = [x y z ones(nr,1)]; % add column of 1's
txyz = (transMat * xyz')';

% Extract for output:
tx = txyz(:,1);
ty = txyz(:,2);
tz = txyz(:,3);


