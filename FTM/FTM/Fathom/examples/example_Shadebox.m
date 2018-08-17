% Example of using the shadebox function
% 
% by David L. Jones, Jun-2012
% 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/shadeBox.mat'

% Variables:
% region = regions of interest to shade
% speed  = velocity in m/sec
% time   = Julian day 1991

% Load the data:
load shadeBox.mat

% Create plot:
plot(time,speed,'b-');

% Label axes:
xlabel('Date (1991)');
ylabel('Speed (m/sec)');

% Shade regions: 
f_shadeBox(region);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Optionally customize the plot:                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List the Gregorian date range:
f_gregorian([min(time) max(time)],1991)
% ans = 
% 
%     '24-May-1991'
%     '01-Jul-1991'

% Get corresponding Julian dates for these Gregorian dates:
% 24-May-1991
% 01-Jun-1991
% 15-Jun-1991
% 01-July-1991
% 
x = [1991 05 24
1991 06 01
1991 06 15
1991 07 01];
% 
jdate = f_julian(x,0,0,1991)
% 
% jdate =
% 
%    144
%    152
%    166
%    182

% Set custom tick mark locations and labels:
txt = {'24-May' '01-Jun' '15-Jun' '01-Jul'};
set(gca,'Xtick',jdate,'XTickLabel',txt)
