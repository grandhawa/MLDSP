% Example of working with NSIDC Sea Ice Coverage data from the South Pole:
% 
% by David L. Jones, July-2012.
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% updated Aug-2012

% File: '.../examples/nt_20110101_f17_nrt_s.bin'

% These data are National Snow & Ice Data Center (NSIDC) sea ice coverage
% data from the south pole downloaded from:
% http://nsidc.org/data/docs/daac/nsidc0081_ssmi_nrt_seaice.gd.html

% Import NSIDC binary file and plot Sea Ice coverage percentages:
A = f_NSIDCimport('nt_20110101_f17_nrt_s.bin');

% Specify site coordinates:
lonLat =[ -26.857 -72.36
   -40.718 -72.243
   -54.24 -69.993
   -34.315 -64.677
   -28.409 -61.516
   -83.328 -72.137
   -91.417 -71.066
   -81.629 -67.571
   -74.404 -64.566
   -146.13 -72.337];

% Plot sites on Sea Ice map:
m_line(lonLat(:,1),lonLat(:,2),'marker','O','markersize',6,'MarkerFaceColor',...
   'w','MarkerEdgeColor','none','LineStyle','none');

% Interpolate Sea Ice values at specific coordinates:
seaIce = f_NSIDCinterp(A,lonLat(:,1),lonLat(:,2))

% seaIce =
% 
%           100
%           100
%           100
%        21.235
%             0
%           100
%           100
%             0
%             0
%        48.354
%        
