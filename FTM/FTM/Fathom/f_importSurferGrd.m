function [lon,lat,z] = f_importSurferGrd(fname,pflag);
% - import Surfer GRID (*.grd) file
%
% USAGE: [lon,lat,z] = f_importSurferGrd('fname',pflag);
%
% fname = name of Surfer *.grd blanking file
% pflag = diagnostic plot flag (default = 0)
%
% lon = longitude or x coordinates
% lat = latitude  or y coordinates
% z   = depth or other
%
% SEE ALSO f_importSurferBln

% -----Notes:-----
% This program is used to import a Golden Software's 'Surfer for 
% Windows' GRID (*.grd) file into Matlab in a format that can be used by
% the M_Map Toolbox. The Surfer Grid file must be exported from Surfer as
% an ASCII XYZ file.

% -----Author(s): -----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input and import file:-----
if (nargin < 2), pflag = 0; end; % no plot by default

if (exist(fname,'file')==0)
   error(['File ' fname ' not found! Check path or filename.']);
end

% Read in data file:
rawData = textread(fname,'','delimiter',' ','headerlines',0);

if size(rawData,2) ~= 3
   error('Incorrect file format for input...check # of columns!');
end
% ----------------------

% -----Parse variables & Reshape matrix:-----
lon = unique(rawData(:,1));
lat = unique(rawData(:,2));
z   = rawData(:,3);
clear rawData; % free up memory
z = reshape(z,size(lon,1),size(lat,1));
z = z';

% Make Diagnostic plot:
if (pflag>0)
   pcolor(lon,lat,z);
   axis equal;
   shading interp;
end;
