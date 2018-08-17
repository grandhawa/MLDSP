function [lon,lat,sst,cmap,h] = f_plotUSGS
% - plot USGS Coastwatch SST image with M_Map Toolbox
%
% USAGE: [lon,lat,sst,cmap,h] = f_plotUSGS;
%
% -----Output:-----
% lon =  vector of longitudes 
% lat =  vector of latitudes
% sst =  indexed-color values representing
%        Sea Surface Temperatures (0-256)
% cmap = image colormap
% h    = graphics handle to surf object

% -----Notes:-----
% A) This function plots a USGS SST 'tif' image using the M_Map Toolbox.
% It assumes the imported file is in Mercator projection, covers an
% area from 22.00°N to 31.61°N and from 79.000°W to 92.658°W, and is
% 960 x 1216 pixels.
%
% B) Matlab 'patches' are created via the SURF function, so this method
% is SLOWER than simply reading/viewing a bitmap file, but allows for much
% more control over the final map.
%
% C) Normally this function should be called with NO output arguments, e.g.,
% >> f_plotUSGS;
%
% D) Edit the lat/long extents in the M_PROJ and M_GRID commands below to crop
% the map.
%
% E) You can easily add additional M_Map commands before the HOLD OFF command
% below to overlay symbols, lines, etc.
%
% Example USGS SST imagery is available from
% http://catbert.er.usgs.gov/east_gulf/html/sst.html

% -----Dependencies:-----
% m_2D_surf.m by Brian Farrelly<brian.Farrelly@nho.hydro.com>
% M_Map Toolbox

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2003: fixed bugs in file import, made compatible with R13

% Import image:
[fileVar,pathVar] = uigetfile('*.tif', 'Specify file to import');
%nameVar = fileVar(1,1:(size(fileVar,2) - 4)); % strip extension off for version 5
nameVar = [pathVar fileVar];
[sst,cmap] = imread(nameVar,'tif');

% Setup x & y axes:
lat = linspace(22.00,31.61,960);
lon = linspace(-92.658,-79.000,1216);

% Format array to be compatible with Matlab's surf:
sst = flipud(double(sst));

% Setup up projection:
m_proj('mercator','longitudes',[-92 -79],'latitudes',[22 31]);

hold on;

% Plot SST:
m_2D_surf(lon,lat,sst); % call Farrelly's SURF routine
colormap(cmap);         % use colormap derived from imported 'tif' file
m_grid('box','fancy','fontsize',8,'linestyle','none','xtick',[-92:2:-79],'ytick',[22:2:31]);

% add additional M_Map commands here before issuing the 'hold off' command;
hold off;

% to save this figure as a, say, 355 dpi 'tif' file:
% print -dtiff -r355 'filename.tif'
