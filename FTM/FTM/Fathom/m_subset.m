function m_subset(fnameIN,fnameOUT,lon,lat)
% - create a subset of an M_Map usercoast file
%
% USAGE: m_subset('fnameIN','fnameOUT',lon,lat);
%
% fnameIN  = name of source file to subset
% fnameOUT = name of destination file
% lon      = min/max of longitude (e.g., lon = [-83 -82]);
% lat      = min/max of latitude  (e.g., lat = [27.25 28.4]);
%
% SEE ALSO: f_importSurfer

% -----References:-----
% Portions of this code are from the comments in mu_coast.m
% from Rich Pawlowicz's <rich@ocgy.ubc.ca> M_Map Toolbox
% available from:
% http://www2.ocgy.ubc.ca/~rich/map.html

% -----Author:-----
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Check input:-----
% Don't overwrite source file:
if isequal(fnameIN,fnameOUT)
   error('Source and destination files are the same!')
end

% Don't overwrite existing file:
if exist(fnameOUT,'file')
   error('Destination file already exists!')
end
% ----------------------

% Load source file:
eval(['load ' fnameIN ' ncst']);

% Define coordinates of rectangle specifying subset (clockwise):
subX = [lon(1) lon(1) lon(2) lon(2) lon(1)]';
subY = [lat(1) lat(2) lat(2) lat(1) lat(1)]';

fprintf('\nSearching for intersections of polygons...\n')
[x,y] = polybool('intersection',ncst(:,1),ncst(:,2),subX,subY);
clear ncst;

fprintf('\nCalculating areas of polygons...\n')
ncst = [[NaN NaN];[x y];[NaN NaN]]; % begin/terminate with NaN's

% ----- This section from mu_coast.m: -----
k=[find(isnan(ncst(:,1)))]; % get indices of NaN's

Area = zeros(length(k)-1,1);
for i = 1:length(k)-1,
   x  = ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
   y  = ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
   nl = length(x);
   Area(i) = sum( diff(x).*(y(1:nl-1)+y(2:nl))/2 );
end

% Area should be >0 for land, and <0 for lakes and inland seas.
% -----------------------------------------

% Sort descending by Area:
fprintf('\nSorting polygons by area...\n');
[Area,ncst,k] = f_sortMap(Area,ncst,k);

% Save subset to destination file:
eval(['save ' fnameOUT ' ncst k Area']);
fprintf(['\nSubset saved as ' fnameOUT '\n'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           SUBFUNCTION                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area2,ncst2,k2] = f_sortMap(Area,ncst,k)
% USAGE [Area2,ncst2,k2] = f_sortMap(Area,ncst,k);

% This function is used to sort map data imported with f_importSurfer.m
% by Area according to the suggestions in mu_coast in the M_Map toolbox

% Remove trailing NaN's:
k = k(1:(end-1));
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
