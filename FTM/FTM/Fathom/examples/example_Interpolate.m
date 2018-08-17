% Example of interpolating values for oceanographic data
% 
% by David L. Jones, July-2012.
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/sst.dat'

% These data are sea surface temperature (SST) values for corresponding
% lat/lon coordinate pairs from the western Antartic Peninsula saved in an
% ASCII spaced-delimited text file. These data are courtsey Paul Suprenand.

% Load data:
raw = load('sst.dat');

% Parse data:
sst = raw(:,1);
lat = raw(:,2)*-1; % change sign of latitudes
lon = raw(:,3)*-1; % change sign of longitudes 
clear raw;

% Plot coordinates to visualize geographic range:
plot(lon,lat,'b.');

% Construct the interpolant:
F = TriScatteredInterp(lon,lat,sst);

% Examine the 100th value from dataset:
[lon(100) lat(100) sst(100)]
% ans =
% 
%        -88.187       -57.313        6.482

% Test interpolate by querying a location with a known SST value coordinates:
qSST = F(-88.187,-57.313)
% 
% qSST = 6.482

% Interpolate SST value at a location where one doesn't already exist:
qSST = F(-70.533,-64.577)
% 
% qSST = 0.416
