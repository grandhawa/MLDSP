% Example of a Vector Plot of Wind Stress
% 
% by David L. Jones, 2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/vecPlot.mat'

% This f_vecPlot function is used to plot time series of wind or current
% meter velocity vectors using Matlab?s quiver function. This function is
% necessary in order to obtain vectors that have the proper length and
% angle of rotation. An optional scaling factor can be applied allowing the
% user control over the amount of overlap among vectors and/or the scaling
% of vectors relative to the overall time series. The X-axis is scaled
% accordingly. The Y-axis allows easy, visual interpretation of vector
% length.
% 
% U,V components of velocity vectors can be extracted from data
% specifying only Speed and Direction using f vecUV.

% Load the data:
load vecPlot.mat
% 
% Create diagram:
figure;
set(gcf,'color','w'); % set bg color to white
f_vecPlot(date,u,v,20,'tau',[140 182]);
f_shadeBox(subsets,20);
