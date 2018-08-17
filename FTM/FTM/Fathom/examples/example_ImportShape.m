% Example of importing an ArcView shape file
% 
% by David L. Jones, Mar-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/strata.shp'

% Import shape file:
[ncst,k,Area] = f_importShapefile('strata.shp');

% Save as an M_Map usercoast file:
save strata.mat ncst k Area;
clear ncst k Area;

% -----Create a Map:-----
figure;               % open new figure window
set(gcf,'color','w'); % set bg color so printed/exported lakes will be white

% Create a base map:
m_proj('mercator','longitudes',[-83 -82],'latitudes',[27.25 28.4]);
m_usercoast('tampaBay','patch',[0.5 0.80 0.5],'edgecolor','none');

% Plot strata:
% note: bottom types must be in work path!
m_usercoast('strata.mat','patch',[1 1 1]*0.85,'edgecolor','none');

% Adjust stacking order, so smaller patches aren't obscured:
h = get(gca,'Children');  % get stacking order of handles
[null,h2] = f_figArea(h); % sort by area
set(gca,'Children',h2);   % re-stack

% Delete 'holes' (i.e., white patches)
handle = findobj(gca,'Type','patch','FaceColor','w');
delete(handle);

% Complete Base Map:
m_grid('box','fancy','fontsize',8,'linestyle','none','xtick',[-83:15/60:-82],...
   'ytick',[27+(15/60):15/60:28+(15/60)]);
% -----------------------
