function [ncst,k,Area] = f_importSurferBln(fname);
% - import Surfer *.bln blanking file for M_Map
%
% USAGE: [ncst,k,Area] = f_importSurferBln('fname');
%
% fname = name of Surfer *.bln blanking file
%
% ncst,k,Area = required components for M_Map usercoast.mat file
%
% SEE ALSO: f_importSurferGrd

% -----Notes:-----
% This program is used to import a Golden Software's 'Surfer for 
% Windows' *.bln file into Matlab in a format that can be used by
% the M_Map Toolbox. Since Surfer can import ArcView Shape files,
% this is a good way to get *.shp data into Matlab. The
% suggested procedure is to turn off the display of axes in Surfer,
% then export the data as a *.bln file (select the option "Break Apart
% Compound Areas"), import using f_importSurfer,
% save as a "usercoast" file, and plot using m_usercoast.m.
%
% Example code to plot as filled patches using M_Map toolbox:
% m_proj('mercator','longitudes',[-81 -80],'latitudes',[24.5 25.5]);
% m_usercoast('fmri.mat','patch',[0 0 0],'edgecolor','none');
% m_grid('box','fancy','fontsize',8,'linestyle','none','xtick',[-81:-80],'ytick',[24.5:25.5]);

% -----Details:-----
% Be sure to save variables ncst,k, and Area as a 'usercoast.mat' file

% -----References:-----
% Portions of this code are from the comments in mu_coast.m
% from Rich Pawlowicz's <rich@ocgy.ubc.ca> M_Map Toolbox
% available from:
% http://www2.ocgy.ubc.ca/~rich/map.html

% -----Author(s): -----
% by Dave Jones, Aug-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2003: added NaN's on last line to be consistent with M_Map format
%           added subfunction f_sortMap, import file
% Mar-2003: renamed as f_importSurferBln

% -----Check input and import file:-----
if (exist(fname,'file')==0)
   error(['File ' fname ' not found! Check path or filename.']);
end

% Read in data file:
bln = textread(fname,'','delimiter',',','headerlines',0);

if size(bln,2) ~= 2
   error('Incorrect file format for input...check # of columns!');
end
% ----------------------

ind = find(bln(:,2) == 0); % get indices of 0's marking beginning of each polygon
bln(ind,1) = NaN;          % replace with NaN's separating closed polygons
bln(ind,2) = NaN;

ncst = [bln;[NaN NaN]]; % terminate with NaN's
clear bln;

% ----- This section from mu_coast.m: -----
k=[find(isnan(ncst(:,1)))]; % get indices of NaN's

Area = zeros(length(k)-1,1);
for i = 1:length(k)-1,
   x  = ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
   y  = ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
   nl = length(x);
   Area(i) = sum( diff(x).*(y(1:nl-1)+y(2:nl))/2 );
end;

% Area should be >0 for land, and <0 for lakes and inland seas.
% -----------------------------------------

% Sort descending by Area:
fprintf('\nSorting objects by Area...please wait \n\n');
[Area,ncst,k] = f_sortMap(Area,ncst,k);

fprintf('Be sure to save ncst,k,Area as a usercoast.mat file for later use \n\n');


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
