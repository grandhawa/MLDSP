function jdate = f_julian(x,zone,daylight,yr)
% - convert date vector to Julian date
%
% USAGE: jdate = f_julian(x,zone,daylight,yr)
%
% x        = input matrix specifying dates     [yyyy mm dd h m s]
% zone     = time zone                         (default = 0)
% daylight = correct for Daylight Savings Time (default = 0)
% yr       = reference year                    (defaults to earliest in x)
% 
% jdate    = Julian date
%
% See also: f_gregorian, f_isDST, f_leapYear

% -----Notes:-----
% This program can span multiple years, with dates after
% the first year > 365 (or 366).
%
% The time zone for US East Coast = -5.
%
% This program calls f_isDST to correct for Daylight Savings Time,
% which currently only supports the US definition.

% ----- Author(s): -----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2009: added support for user-supplied reference year
% Jun-2012: updated documentation

% -----Check input and set defaults:-----
if (nargin < 2), zone     = 0; end; % default no correction for time zone
if (nargin < 3), daylight = 0; end; % default no correction for DST
if (nargin < 4), yr       = 0; end; % default no correction for DST

if sum(size(yr),2)~=2, error('YR must be a scalar!'); end

% Pad input with 0's if < 6 columns:
[nr,nc] = size(x);
if (nc < 6)
	x(nr,nc+1:6) = 0;
end
% ---------------------------------------

% Convert to serial date (days since 01-Jan-0000):
serialDate = datenum([x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6)]);

% Apply time zone correction:
if (zone ~= 0)
	serialDate = serialDate + (zone*(1/24));
end

% add 1 hr to dates during Daylight Savings:
if (daylight>0)
	[y,m,d,h] = datevec(serialDate);
	dst       = f_isDST([y m d h]);
	serialDate(find(dst==1)) = serialDate(find(dst==1)) + (1/24);
end

% Get starting date (31-Dec of previous year):
if (yr>0)
   startDay = datenum([yr-1,12,31,0,0,0]);          % use reference year provided
else
   startDay = datenum([min(x(:,1))-1,12,31,0,0,0]); % input dates

end

% Rescale (days since 31-Dec of prevous year):
jdate = (serialDate - startDay);

