% Examples of creating various maps using the M_Map Toolbox:
% 
% by David L. Jones, July-2012.
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Mollweide Equal-Area Projection Map of the western Antartic Peninsula:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% File: '.../examples/sst.dat'

% M_Map setup:
figure;
set(gcf,'color','w'); % set bg color so printed/exported lakes will be white

% Create a base map:
m_proj('mollweide','longitudes',[-79 -54],'latitudes',[-73 -61]);
% m_coast('patch',[1 1 1]*0.75,'edgecolor','none');
% m_gshhs_h('save','WAP_h.mat'); % save region of interest
m_usercoast('WAP_h.mat','patch',[1 1 1]*0.75,'edgecolor','none');

% -----Complete Base Map:-----
m_grid('box','fancy','fontsize',8,'linestyle',':',...
   'xtick',(-78:4:-54),'ytick',(-72:4:-60));

% Add scale bar:
m_scaleBar(-78.25,-62,250,4,8)
