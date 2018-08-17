function [] = f_cpsParse(raw)
% - GUI to parse LA-ICP-MS transient signal data into signal/background subsets
%
% USAGE: f_cpsParse(raw)
%
% INPUT:
% raw = structure of raw signal data with the following fields:
%  .s     = seconds
%  .cps   = signal intensity (counts per second)
%  .txt   = cell array of element labels for columns of 'cps'
%  .iso   = corresponding isotope number
%  .hh    = hour of acquisition
%  .mm    = minute of acquisition
%  .gDate = date of acquition
%
% OUTPUT:
% Subset of 'raw' input structure with the following fields:
%  .s     = seconds
%  .cps   = signal intensity (counts per second)
%  .txt   = cell array of element labels for columns of 'cps'
%  .iso   = corresponding isotope number
%  .hh    = hour of acquisition
%  .mm    = minute of acquisition
%  .gDate = cell array of Gregorian date of acquition
%  .idx   = index to corresponding elements of input 'raw' structure
%  .bg    = structure with the corresponding background data
%
% SEE ALSO: f_importXL, f_cps2ppm, f_cpsPlot

% -----Notes:-----
% This function is used to parse a time series of laser ablation transient
% signal data into separate components representing the signal/background
% portions of individual ablation samples. The input data (RAW) must be a
% structure created from the 'f_importXL' function, which imports data from a
% Perkin-Elmer Elan 'XL' format data file. The GUI created by this function
% displays a plot of the time series data that allows the user to parse by
% selecting specific regions. Use the 'Pan' or 'Zoom' tools to navigate and
% focus on specific portons of the time series. Use the 'Select' tool to
% highlight a region within the time series that corresponds with the SIGNAL
% portion of a sample, add an appropriate name (e.g., spot_01) within the GUI's
% 'Variable Name' box, then click the 'Signal' button to save this portion of
% the data in a separate variable in the Matlab workspace. Next, use the
% 'Select' tool to highlight the region of the time series that corresponds with
% the BACKGROUND portion of the sample you just exported, and click the
% 'Background' button. This adds a *.bg field to the signal variable you
% just created (i.e., to the variable named in the 'Variable Name' box).
% Once you've parse the entire data set, use the 'f_cps2ppm' function to
% process the data.
% 
% There appears to be a slight bug in the way Matlab creates the GUI, so
% when the window is first drawn and the user selectes the ZOOM tool for
% the fist time, the region of interest (ROI) drawn by the ZOOM tool is
% invisible, though the ZOOM tool still works. To avoid this, simply use
% the 'Select' tool to draw an arbitraty ROI before using the ZOOM tool for
% the first time.
% 
% 
% If your screen is small, you may need to edit the size of the GUI window.
% This works on a 14" laptop:
% S.fh = figure('units','pixels',...
%    'position',[75 100 1024 600],...
%    'menubar','none',...
%    'numbertitle','off',...
%    'name',nameVar,...
%    'resize','off');
% 
% % Set axes:
% S.ax = axes('units','pixels',...
%    'position',[40 160 900 660],...
%    'fontsize',8);

% -----Author:-----
% by David L. Jones, Aug-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2010: added support for time stamp, updated documentation.
% Jun-2011: added Background button; add date of acquisition
% Apr-2012: edited documentation

% -----Check input & check defaults:-----
% Extract data:
s       = raw.s;
cps     = raw.cps;
txt     = raw.txt;
iso     = raw.iso;
hh      = raw.hh;
mm      = raw.mm;
gDate   = raw.gDate;
clear raw;

% Get time interval
intTime = s(2) - s(1);
nameVar = sprintf('Select region of ''%s'' (interval = %g sec):',inputname(1), intTime);
% ---------------------------------------

% -----Setup GUI:-----
S.fh = figure('units','pixels',...
   'position',[75 100 1024 850],...
   'menubar','none',...
   'numbertitle','off',...
   'name',nameVar,...
   'resize','off');

% Set axes:
S.ax = axes('units','pixels',...
   'position',[40 160 900 660],...
   'fontsize',8);

% Initialize structure with data:
S.s       = s;       % dependent timeseries variable
S.cps     = cps;     % independent timeseries variable
S.txt     = txt;     % cell array of element text labels
S.iso     = iso;     % corresponding isotope number
S.hh      = hh;      % hour of acquisition
S.mm      = mm;      % min  of acquisition
S.gDate   = gDate;   % date of acquisition
S.xlm     = [];      % x-axis limit
S.ylm     = [];      % y-axis limit
S.hl      = [];      % handle to plot of total counts
S.hl2     = [];      % handle to plot of individual counts
S.hp      = [];      % handle to patch
S.min     = [];      % initialize
S.max     = [];

% Plot plot total counts:
S.hl = plot(S.s,sum(S.cps,2));
set(S.ax,'Yscale','log','YLimMode','auto');
xlabel('Seconds')
ylabel('Signal Intensity (cps)');

% Prevent lines from responding to mouse clicks:
set(S.hl,'Hittest','off');

% Update structure:
S.xlm = get(S.ax,'xlim'); % x-axis limit
S.ylm = get(S.ax,'ylim'); % y-axis limit

% 'Min' text label:
S.tx(1) = uicontrol('style','tex',...
   'unit','pix',...
   'posit',[50 70 120 25],...
   'backg',get(S.fh,'color'),...
   'fontsize',12,'fontweight','bold',...
   'string','Min:');

% 'Max' text label:
S.tx(2) = uicontrol('style','tex',...
   'unit','pix',...
   'posit',[190 70 120 25],...
   'backg',get(S.fh,'color'),...
   'fontsize',12,'fontweight','bold',...
   'string','Max:');

% 'Variable Name' text label:
S.tx(3) = uicontrol('style','tex',...
   'unit','pix',...
   'posit',[330 70 240 25],...
   'backg',get(S.fh,'color'),...
   'fontsize',12,'fontweight','bold',...
   'string','Variable Name:');

% 'Min' edit box::
S.ed(1) = uicontrol('style','edit',...
   'unit','pix',...
   'position',[50 50 120 25],...
   'backg','w',...
   'fontsize',10,'fontweight','bold' );

% 'Max' edit box:
S.ed(2) = uicontrol('style','edit',...
   'unit','pix',...
   'position',[190 50 120 25],...
   'backg','w',...
   'fontsize',10,'fontweight','bold');

% Variable name ('idxName') edit box:
S.ed(3) = uicontrol('style','edit',...
   'unit','pix',...
   'position',[330 50 240 25],...
   'fontsize',10,...
   'backg','w',...
   'fontsize',10,'fontweight','bold',...
   'string','idxName');

% 'Signal' push button:
S.pb(1) = uicontrol('style','push',...
   'unit','pix',...
   'position',[590 50 120 25],...
   'fontsize',10,'fontweight','bold',...
   'string','Signal');

% 'Total CPS' push button:
S.pb(2) = uicontrol('style','push',...
   'unit','pix',...
   'position',[720 50 120 25],...
   'fontsize',10,'fontweight','bold',...
   'string','Total CPS');

% 'Background' push button:
S.pb(3) = uicontrol('style','push',...
   'unit','pix',...
   'position',[590 25 120 25],...
   'fontsize',10,'fontweight','bold',...
   'string','Background');


%-----Radio Button Group:-----
% Note position of radio buttons are relative to button group
S.bg(1) = uibuttongroup('units','pix',...
   'pos',[850 35 90 100]);
S.rd(1) = uicontrol(S.bg(1),'style','rad',...
   'unit','pix',...
   'position',[10 65 70 30],...
   'fontsize',10,'fontweight','bold',...
   'string','Select');
S.rd(2) = uicontrol(S.bg(1),'style','rad',...
   'unit','pix',...
   'position',[10 35 70 30],...
   'fontsize',10,'fontweight','bold',...
   'string','Zoom');
S.rd(3) = uicontrol(S.bg(1),'style','rad',...
   'unit','pix',...
   'position',[10 5 70 30],...
   'fontsize',10,'fontweight','bold',...
   'string','Pan');
% -----------------------

% Set callback functions:
set(S.ax,'buttondownfcn',{@ax_bdfcn,S});      % axes
set(S.ed(1),'callback',{@ed_call,S});         % Min field
set(S.ed(2),'callback',{@ed_call,S});         % Max field
set(S.pb(1),'callback',{@pb_call,S});         % Signal button
set(S.pb(2),'callback',{@pb2_call,S});        % Total CPS button
set(S.pb(3),'callback',{@pb3_call,S});        % Background button
set(S.bg(1),'SelectionChangeFcn',{@bg_call}); % Radio button group


% -----Buttondownfcn for axes:-----
function [] = ax_bdfcn(varargin)
[h,S] = varargin{[1,3]};  % extract the calling handle and structure.

delete(S.hp); % delete any previous patch

pnt1      = get(h,'CurrentPoint'); % button down detected
finalRect = rbbox;                 % rubberban box returns figure units
pnt2      = get(h,'CurrentPoint'); % button up detected
pnt1      = pnt1(1,1:2);           % extract x and y
pnt2      = pnt2(1,1:2);
p1        = min(pnt1,pnt2);        % calculate locations
offset    = abs(pnt1-pnt2);        % and dimensions
x         = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
% y       = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
y         = [S.ylm(1) S.ylm(1) S.ylm(2) S.ylm(2) S.ylm(1)]; % extend to Y-limits

% Replace X with nearest observed value of s:
for i = 1:2
   xDiff    = abs(S.s-x(i)); % difference between actual time and selected value;
   idx      = find(xDiff == min(xDiff));
   S.idx(i) = max(idx);      % for ties, select the larger index
   x(i)     = S.s(S.idx(i));
end
x(3:5) = [x(2) x(1) x(1)];

% Plot selected region as semi-transparent patch:
S.hp = patch(x,y,'k');
% Semi-transparent, no mouse response:
set(S.hp,'EdgeColor','none','FaceColor',[1 1 1]*0.80,'FaceAlpha',0.25,'Hittest','off');

% Save bounds of selected time range:
S.min = min(x);
S.max = max(x);

% Update edit boxes:
set(S.ed(1),'str',num2str(S.min,'%g')); % x min
set(S.ed(2),'str',num2str(S.max,'%g')); % x max

% Make updated S available to other callbacks:
set(S.ax,'buttondownfcn',{@ax_bdfcn,S}); % axes
set(S.ed(1),'callback',{@ed_call,S});    % Min field
set(S.ed(2),'callback',{@ed_call,S});    % Max field
set(S.pb(1),'callback',{@pb_call,S});    % Export button
set(S.pb(2),'callback',{@pb2_call,S});   % Total CPS button
set(S.pb(3),'callback',{@pb3_call,S});   % Background button
% ---------------------------------


% -----Callback for Min/Max edit boxes:-----
function [] = ed_call(varargin)
S = varargin{numel(varargin)}; % extract structure

% Get values from edit box:
x(1) = str2double(get(S.ed(1), 'String')); % xMin
x(2) = str2double(get(S.ed(2), 'String')); % xMax

% Check range of input:
if (x(2) <= x(1)), error('X Max should be > X Min!'); end

% Replace x with nearest observed value of s:
for i = 1:2
   xDiff    = abs(S.s-x(i)); % difference between actual time and selected value;
   idx      = find(xDiff == min(xDiff));
   S.idx(i) = max(idx);      % for ties, select the larger index
   x(i)     = S.s(S.idx(i));
end

% Update structure:
S.min = x(1);
S.max = x(2);
set(S.hp,'Xdata',[S.min S.max S.max S.min S.min]); % X-boundaries of patch object

% Update edit boxes:
set(S.ed(1),'str',num2str(S.min,'%g')); % x min
set(S.ed(2),'str',num2str(S.max,'%g')); % x max

% Make updated S available to other callbacks:
set(S.ed(1),'callback',{@ed_call,S});    % Min field
set(S.ed(2),'callback',{@ed_call,S});    % Max field
set(S.pb(1),'callback',{@pb_call,S});    % Export button
set(S.pb(2),'callback',{@pb2_call,S});   % Total CPS button
set(S.pb(3),'callback',{@pb3_call,S});   % Background button


% -----SelectionChangeFnc for Radio Buttons:-----
function [] = bg_call(source,event)
btn = (get(event.NewValue,'String'));
switch btn
   case 'Select'
      pan off;
      zoom off;
   case 'Zoom'
      pan off;
      zoom on;
   case 'Pan'
      zoom off;
      pan on;
end


% -----Callback for 'Signal' push button:-----
function [] = pb_call(varargin)
S     = varargin{numel(varargin)}; % extract structure
vName = get(S.ed(3), 'String');    % get name of variable to export to

if ~isempty(vName)
   % 'assignin' can't work with structures, so:
   assignin('base','assignin_temp_s',S.s(S.idx(1):S.idx(2)));
   assignin('base','assignin_temp_cps',S.cps(S.idx(1):S.idx(2),:));
   assignin('base','assignin_temp_txt',S.txt);
   assignin('base','assignin_temp_iso',S.iso);
   assignin('base','assignin_temp_hh',S.hh);
   assignin('base','assignin_temp_mm',S.mm);
   assignin('base','assignin_temp_gDate',S.gDate);
   assignin('base','assignin_temp_idx',S.idx); % keep a record of indices
   evalin('base',[vName '.s       = assignin_temp_s;'])
   evalin('base',[vName '.cps     = assignin_temp_cps;'])
   evalin('base',[vName '.txt     = assignin_temp_txt;'])
   evalin('base',[vName '.iso     = assignin_temp_iso;'])
   evalin('base',[vName '.hh      = assignin_temp_hh;'])
   evalin('base',[vName '.mm      = assignin_temp_mm;'])
   evalin('base',[vName '.gDate   = assignin_temp_gDate;'])
   evalin('base',[vName '.idx     = assignin_temp_idx;'])
   evalin('base','clear assignin_temp*')
end


% -----Callback for 'Background' push button:-----
function [] = pb3_call(varargin)
S     = varargin{numel(varargin)}; % extract structure
vName = get(S.ed(3), 'String');    % get name of variable to export to

if ~isempty(vName)
   % 'assignin' can't work with structures, so:
   assignin('base','assignin_temp_s',S.s(S.idx(1):S.idx(2)));
   assignin('base','assignin_temp_cps',S.cps(S.idx(1):S.idx(2),:));
   assignin('base','assignin_temp_txt',S.txt);
   assignin('base','assignin_temp_iso',S.iso);
   assignin('base','assignin_temp_hh',S.hh);
   assignin('base','assignin_temp_mm',S.mm);
   assignin('base','assignin_temp_gDate',S.gDate);
   assignin('base','assignin_temp_idx',S.idx); % keep a record of indices
   evalin('base',[vName '.bg.s       = assignin_temp_s;'])
   evalin('base',[vName '.bg.cps     = assignin_temp_cps;'])
   evalin('base',[vName '.bg.txt     = assignin_temp_txt;'])
   evalin('base',[vName '.bg.iso     = assignin_temp_iso;'])
   evalin('base',[vName '.bg.hh      = assignin_temp_hh;'])
   evalin('base',[vName '.bg.mm      = assignin_temp_mm;'])
   evalin('base',[vName '.bg.gDate   = assignin_temp_gDate;'])
   evalin('base',[vName '.bg.idx     = assignin_temp_idx;'])
   evalin('base','clear assignin_temp*')
end


% -----Callback for 'Total CPS' push button:-----
function [] = pb2_call(varargin)
S = varargin{numel(varargin)}; % extract structure

if isempty(S.hl2)
   % Plot plot individual counts:
   hold on;
   S.hl2 = plot(S.s,S.cps);
   hold off;
   
   % Prevent lines from responding to mouse clicks:
   set(S.hl2,'Hittest','off');
else
   delete(S.hl2);
   S.hl2 = [];
end


% Update Y-values of patch to match new axis limits:
if ~isempty(S.hp)
   S.ylm = get(S.ax,'ylim'); % current y-axis limit
   y     = [S.ylm(1) S.ylm(1) S.ylm(2) S.ylm(2) S.ylm(1)]; % extend to Y-limits
   set(S.hp,'YData',y);
end

% Make updated S available to other callbacks:
set(S.ax,'buttondownfcn',{@ax_bdfcn,S}); % axes
set(S.ed(1),'callback',{@ed_call,S});    % Min field
set(S.ed(2),'callback',{@ed_call,S});    % Max field
set(S.pb(1),'callback',{@pb_call,S});    % Signal button
set(S.pb(2),'callback',{@pb2_call,S});   % Total CPS button
set(S.pb(3),'callback',{@pb3_call,S});   % Background button

