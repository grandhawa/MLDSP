function [jdate,u,v,ws,wd] = f_windCMAN(fname,zone,daylight)
% - process CMAN (or NDBC) historical wind data
%
% USAGE: [jdate,u,v,ws,wd] = f_windCMAN('fname',zone,daylight)
%
% fname    = name of CMAN historical data file
% zone     = time zone                         (default = 0)
% daylight = correct for Daylight Savings Time (default = 0)
%
% jdate = Julian date (optionally corrected for time zone & DST)
% u     = East wind component 
% v     = North wind component
% ws    = wind speed (m/s)
% wd    = direction wind is blowing to (North = 0, East = 90 degrees)
%
% SEE ALSO: f_windstress, f_vecPlot, f_ekmanDepth, f_vecRot

% ----- Notes: -----
% This function is used to process wind data from a CMAN historical
% data file downloaded from, for example:
% http://www.ndbc.noaa.gov/station_history.phtml?station=smkf1
%
% The format of these files is as follows, with 999 & 99.0 indicating
% missing values:
% YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS
% 92 01 01 00 353 06.2 06.6 99.00 99.00 99.00 999 1019.2  20.0  24.1 999.0 99.0
%
% Wind direction is converted "From" to "To" and a geographic coordinate
% system is used to calculate U,V (= East,North) components of wind. These
% data can then be used in f_vecRot (with the GEO=1 option) to rotate the
% coordinate system clockwise by the compass angle of the local isobaths
% (i.e., the shoreline) to provide the cross-shore (u) and alongshore (v)
% wind components.
%
% Use f_vecPlot to plot the un-rotated u,v data output by this function:
% North = up, East = right, West = left, South = Down; or rotate into
% alongshore/cross-shore components (via f_vecRot) in which case: +V = up
% and +U = right. 

% ----- References: -----
% with help from Cynthia Yeung
% and news://comp.soft-sys.matlab

% -----Author:-----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Dec-2002: added convert 2 digit yr to 4 digit; fixed documentation,
%           mark vs. remove missing data
% Oct-2003: fixed bug in conversion of wind direction and vector rotation,
%           updated documentation
% Oct-2005: removed rotaton by THETA degrees, better to use f_vecRot;
%           removed adding 90 degrees to bring NORTH up as this is done
%           correctly now by f_vecUV with the GEO=1 option now; round julian   
%           dates to 4 decimal places; added to documentation.

% ----- Check input & set defaults: -----
if (nargin < 2), zone     = 0; end; % default no correction for time zone
if (nargin < 3), daylight = 0; end; % default no correction for DST

if (exist(fname,'file')==0)
	error(['File ' fname ' not found! Check path or filename.']);
end
% ---------------------------------------

% read in data file:
x = textread(fname,'','delimiter',' ','headerlines',1);

% Mark missing data:
null_wd = find(x(:,5)==999);  % wind speed missing
null_ws = find(x(:,6)==99.0); % wind direction missing
x(unique([null_wd' null_ws']'),5:6) = NaN;

% Parse data:
yy = x(:,1); % year
mm = x(:,2); % month
dd = x(:,3); % day  in GMT
hh = x(:,4); % hour in GMT
wd = x(:,5); % direction wind is blowing from (North is 0 degrees)
ws = x(:,6); % average wind speed for 2 mins during each hour (m/s)
clear x;

% convert 2 digit years to 4 digits (1980 is pivot year):
yy(find((yy>=80) & (yy<100))) = yy(find((yy>=80) & (yy<100))) + 1900;
yy(find(yy<80)) = yy(find(yy<80)) + 2000;

% Julian date:
jdate = f_julian([yy mm dd hh],zone,daylight);
jdate = f_round(jdate,4); % round to 4 decimal places

% Convert wind direction "From" to "To":
wd = wd + 180;
wd(find(wd>360)) = wd(find(wd>360)) - 360;

% East (u) and North (v) wind components:
[u,v] = f_vecUV(ws,wd,1);
