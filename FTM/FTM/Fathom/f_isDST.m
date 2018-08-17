function dst = f_isDST(date)
% - determine if dates occur during Daylight Savings Time
%
% USAGE: dst = f_isDST(date)
%
% date = matrix of dates having 4 columns [yy mm dd hh]
% dst  = Boolean flag indicating DST (1 = true, 0 = false)
%
% See also: f_julian, f_gregorian, f_leapYear

% -----Notes:-----
% In the US, Daylight Savings Time starts at 2am on the first
% Sunday in April (2am becomes 3am), and ends at 2am on the last
% Sunday in October (2am becomes 1am).
%
% In Europe, Daylight Savings Time starts on the last Sunday
% in March, so slight modification of this code is required.

% -----Author(s):-----
% modified after Matlab code by Randy Poe<rpoe@atl.lmco.com>
% posted to news://comp.soft-sys.matlab (02-July-2002)
%
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (size(date,2)<4)
	error('4 columns specifying [yy mm dd hh] are required')
end

noRows = size(date,1);

% Parse input:
y = date(:,1);
m = date(:,2);
d = date(:,3);
h = date(:,4);

serialDate = datenum([y,m,d,h,zeros(noRows,1),zeros(noRows,1)]);
dow = weekday(serialDate); % day-of-week (Sunday = 1)

dst = zeros(noRows,1); % preallocate results

% Mark dates during DST with 1:
dst(find(m>4   & m<10)) = 1;
dst(find(m==4  & d>7))  = 1;
dst(find(m==4  & d<=7 & dow==1 & h>= 2)) = 1;
dst(find(m==4  & d<=7 & dow>1  & d>dow)) = 1;
dst(find(m==10 & d<25)) = 1;
dst(find(m==10 & d>=25 & dow==1 & h<2)) = 1;
dst(find(m==10 & d>=25 & dow>1 & d-24-dow < 1)) = 1;


