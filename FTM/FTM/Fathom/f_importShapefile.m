function [ncst,k,Area] = f_importShapefile(fname,utm_zone)
% - import ArcView shapefile for M_Map Mapping Toolbox
%
% USAGE: [ncst,k,Area] = f_importShapefile('fname','utm_zone');
%
% fname    = name of ArcView *.shp shape file
% utm_zone = optional UTM zone, if provided UTM converted to LON/LAT
%            (e.g., utm_zone = '17 R')
%
% ncst,k,Area = required components for M_Map usercoast.mat file
%
% SEE ALSO: f_importSurfer, f_utm2deg, m_subset

% -----Notes:-----
% This program is used to import an ArcView shape file into Matlab in a
% format that can be used by the M_Map Toolbox. The resulting variables
% (ncst,K,Area) should be saved as a "usercoast" file, and can be plotted using
% the m_usercoast.m function. Note this function requires both the M_Map
% Mapping Toolbox and Matlab's own Mapping Toolbox.
%
% For converting UTM coordinates to lon/lat this function calls f_utm2deg.m
% which was obtained as UTM2DEG from:
% http://www.mathworks.com/matlabcentral/fileexchange/.
%
% Note this function requires that you specify the UTM zone, which can be
% determine by UTMZONE from the Matlab Mapping Toolbox to determine which
% zone to specify.
%
% See 'exampleImportShape.m' for an example.

% -----Details:-----
% Be sure to save variables ncst,k, and Area as a 'usercoast.mat' file

% -----References:-----
% Portions of this code are from the comments in mu_coast.m
% from Rich Pawlowicz's <rich@ocgy.ubc.ca> M_Map Toolbox
% available from:
% http://www2.ocgy.ubc.ca/~rich/map.html

% -----Author:-----
% by David L. Jones, Mar-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2011: edited to handle changes to f_utm2deg

% -----Check input and import file:-----
if (exist(fname,'file')==0)
   error(['File ' fname ' not found! Check path or filename.']);
end

if (nargin < 2), utm_zone=0; end % default no conversion from UTM to LON/LAT
% ----------------------

% Import shapefile as a 'mapstruct':
S = shaperead(fname);

% Extract lon/lat coordinates:
n    = size(S); % get # elements in mapstruct
ncst = [];      % initialize
for i=1:n
   ncst = [ncst; [S(i).X' S(i).Y']];
end
clear S;

% Optionally convert UTM to Lon/Lat:-----
if ~isequal(0,utm_zone)
   crds = f_utm2deg(ncst(:,1),ncst(:,2),repmat(utm_zone,size(ncst,1),1));
   ncst = [crds.lon crds.lat];
end

% ----- This section from mu_coast.m: -----
k = find(isnan(ncst(:,1))); % get indices of NaN's

Area = zeros(length(k)-1,1);
for i = 1:length(k)-1,
   x       = ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
   y       = ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
   nl      = length(x);
   Area(i) = sum( diff(x).*(y(1:nl-1)+y(2:nl))/2 );
end

% Area should be >0 for land, and <0 for lakes and inland seas.
% -----------------------------------------

% Sort descending by Area:
fprintf('\nSorting objects by Area...please wait \n\n');
[Area,ncst,k] = sub_sortMap(Area,ncst,k);

fprintf('Be sure to save ncst,k,Area as a usercoast.mat file for later use \n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           SUBFUNCTION:                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area2,ncst2,k2] = sub_sortMap(Area,ncst,k)
% USAGE [Area2,ncst2,k2] = sub_sortMap(Area,ncst,k);

% This function is used to sort map data imported with f_importShapefile
% by Area according to the suggestions in mu_coast in the M_Map toolbox

% Remove trailing NaN's:
k    = k(1:(end-1));
ncst = ncst(1:(end-1),:);

% Get indices of blocks of data separated by NaN's:
k2  = k-1;
k2  = [k2(2:end);size(ncst,1)];
idx = [k k2]; % indices specifying [start end] of blocks

% Extract each block separately into a cell array:
for i=1:size(idx,1)
   blocks{i} = ncst(idx(i,1):idx(i,2),:);
end

% Sort blocks by Area, descending:
[Area2,key] = sort(Area);
Area2       = flipud(Area2);
key         = flipud(key);

ncst2 = []; % preallocate
for j=1:size(key,1);
   ncst2 = [ncst2; blocks{key(j)}];
end

% Return trailing NaN:
ncst2 = [ncst2;[NaN NaN]];

% Get indices of NaN's:
k2 = find(isnan(ncst2(:,1))==1);
