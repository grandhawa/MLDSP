function f_gshhs2Shp(lon,lat,fname)
% - export a shapefile from the GSHHS database for use in Surfer/Arcview
%
% USAGE: f_gshhs2Shp([lonMin lonMax],[latMin latMax],'fname.shp');
%
% [lonMin lonMax] = row vector specifying min/max of longitudes
% [latMin latMax] = row vector specifying min/max of latidudes
% 'fname'         = character array specifyin export filename
%
% SEE ALSO: f_importSurfer

% -----References:-----
% After demo in Mapping Toolbox 2006a

% -----Notes:-----
% This script requires the MATLAB Mapping Toolbox in version 2006a.
%
% For some applications, you may not need to process the levels, just
% export the geostructue (e.g., 'geo' below) with shapewrite. Also, this
% code can be edited to optionally allow the user to extract data from the
% other resolutions of gshhs.

% -----Author:-----
% by David L. Jones, Jun-2006
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Create geostruct from the intermediate resolution database:
geo = gshhs('gshhs_i.b',lat,lon);

% See the levels you've imported
levels = [geo.Level];
% Unique(levels) =  [1 2 3];

% Extract coastlines:
level_1 = geo(levels == 1);

% Extract lakes/seas within coastlines
level_2 = geo(levels == 2);

% % Create map to check data:
% figure;
% axesm('mercator', 'MapLatLimit', lat, 'MapLonLimit', lon)
% gridm; mlabel; plabel
% geoshow([level_1.Lat], [level_1.Lon], 'Color', 'blue')
% geoshow([level_2.Lat], [level_2.Lon], 'Color', 'red')
% tightmap

% Copy the level_1 bounding boxes to a 3-D array for convience in looping:
level_1Boxes = reshape([level_1.BoundingBox],[2 2 numel(level_1)]);

for k = 1:numel(level_2)
   for j = 1:numel(level_1)
      % See if bounding boxes intersect:
      if boxesIntersect(level_2(k).BoundingBox, level_1Boxes(:,:,j));
         % See if actual features intersect:
         if ~isempty(polybool('intersection',level_2(k).Lon, level_2(k).Lat,...
               level_1(j).Lon, level_1(j).Lat))
            % Reverse level_2 vertex order before merge to correctly orient
            % inner rings:
            level_1(j).Lon = [level_1(j).Lon fliplr(level_2(k).Lon) NaN];
            level_1(j).Lat = [level_1(j).Lat fliplr(level_2(k).Lat) NaN];
         end
      end
   end
end


% Avoid problem with missing value of 'FormatVersion' field when running
% shapewrite (DLJ):
for i=1:numel(level_1)
   level_1(i).FormatVersion = 3;
end

% Export shapefile (as *.shp, *.shx, and *.dbf):
shapewrite(level_1,fname);

% Validate exported shapefile by importing/plotting:
testShp = shaperead(fname,'UseGeoCoords',true);
figure;
ax = axesm('mercator', 'MapLatLimit', lat, 'MapLonLimit', lon);
set(ax, 'Color', 'cyan')
gridm; mlabel; plabel
geoshow(testShp, 'FaceColor', [0.15 0.8 0.15])
tightmap


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           SUBFUNCTION                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine if bounding boxes intersect (1=true, 0=false):
function res = boxesIntersect(A,B)
res = (~(any(A(2,:) < B(1,:)) || any(B(2,:) < A(1,:))));
