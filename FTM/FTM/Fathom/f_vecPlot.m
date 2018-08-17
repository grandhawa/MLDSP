function f_vecPlot(jdate,u,v,scale,units,jRange)
% - plot time series of velocity vectors
%
% USAGE: f_vecPlot(jdate,u,v,scale,units,jRange)
%
% jdate  = column vector of Julian dates
% u,v    = corresponding vector components
% scale  = scale factor                       (default = 1)
% units  = Y-axis label; e.g., units = 'm/s') (default = none)
% jRange = limits of dates to plot            (default = auto)
%          (e.g., jRange = [min max])
%
% See also: f_julian, f_vecUV, f_shadeBox

% ----- Notes: -----
% This function is used to plot time series of wind or current meter
% velocity vectors using Matlab's QUIVER function. This function is 
% necessary in order to obtain vectors that have the proper length and
% angle of rotation.
%
% An optional scaling factor (SCALE) can be applied allowing the user
% control over the amount of overlap among vectors and/or the scaling of
% vectors relative to the overall time series. The X-axis is scaled
% accordingly. The Y-axis allows easy, visual interpretation of vector
% length. 
%
% U,V components of velocity vectors can be extracted from
% data specifying only Magnitude and Direction using f_vecUV.
%
% MINOR TICKS:
% Display minor tick marks on the x-axis via:
% >> set(gca,'XMinorTick','on');
%
% NORTH ARROWS:
% You may add a North arrow to a plot by using the "Insert Text" tool in
% the Matlab figure window. A nice touch is to use a font from ESRI's
% ArcView or Golden Software's Surfer programs and rotate it as required.
% For example, type a "!", change the font to "GSI North Arrows", increase the
% font size, and rotate as appropriate.
%
% ----- Author: -----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Dec-2003: scaling now operates on the X- vs. Y-axis; Y-axis limits can now be
%           specified 
% Oct-2003: optional angle of rotation (theta), updated documentation
% Nov-2003: updated documentation on 'USAGE'
% Oct-2004: set default line clipping to off
% Jun-2005: removed the following commands to allow use in a subplot: figure;
%           hold on; hold off; set(gcf,'color','w'); 
% Sep-2005: commented out legend command to allow use in a subplot
% Oct-2005: removed rotation angle THETA, better to use f_vecRot (with
%           geo=1)
% Aug-2006: minor edits; also forced setting of XTICK to avoid problems
%           with zooming, scaling, exporting (needed for v.2006a) 

range = 1; % initialize

% ----- Check input & set defaults: -----
if (nargin < 4), scale  = 1;  end; % no scaling by default
if (nargin < 5), units  = []; end; % no units by default
if (nargin < 6), range  = 0;  end; % no range specified

if (scale==0)
   error('You cannot scale vectors by 0');
end

if (size(u,1) ~= size(v,1)) || (size(u,1) ~= size(jdate,1))
   error('U,V, and JDATE must be same size!')
end

set(0,'DefaultLineclipping','off'); % don't clip the data
% ---------------------------------------

nr = size(jdate,1); % # rows

% =========================================================
% figure;
hold on;

% Plot vectors:
quiver(jdate/scale,zeros(nr,1),u,v,0,'.b-');

% Plot base line:
plot([min(jdate/scale) max(jdate/scale)]',[0 0]','k-');

% Adjust aspect ratio for correct angles and lengths:
daspect([1 1 1]);

% Adjust X-axis limits & labels:
if (scale ~= 1)
   
   if (range>0)
      xlim([jRange(1)/scale jRange(2)/scale]); % set x-axis limits
   end
   
   xLabels = get(gca,'xticklabel');  % get tick labels
   xLabels = str2num(xLabels);       % convert to numbers
	set(gca,'XTick',xLabels);         % force setting of tick marks
   xLabels = num2str(xLabels*scale); % recale values
   set(gca,'xticklabel',xLabels);    % replace labels
   
else
   if (range>0)
      xlim([jRange(1)/scale jRange(2)/scale]); % set x-axis limits
   end
end

% Adjust plot appearance:
% set(gcf,'color','w');
set(gca,'TickDir','out');
xlabel('Julian Day');
ylabel(units);
box off;
