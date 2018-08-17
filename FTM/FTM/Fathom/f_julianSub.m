function idx = f_julianSub(jdate,h)
% - get indices to subsample a series of Julian dates
%
% USAGE: idx = f_julianSub(jdate,h);
%
% jdate = column vector of Julian dates
% h     = integer specifying how often to subsample
%         (h = every 6, 12, or 24 hrs)
%
% idx = indices of jdate corresponding to subsample
%
% SEE ALSO f_julian

% -----Notes:-----
% If JDATE is a column vector of Julian dates, with one date every hour,
% suppose you want to subsample a time series every 6 hours. Do the
% following:
%
% > idx     = f_julianSub(jdate,6);
% > sixHour = jdate(idx);

% -----Author:-----
% by Dave Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


jdate = jdate(:); % force column vector

% determine fractional part of hourly julian dates:
switch h
   case 6
      jh = 0.25;
   case 12
      jh = 0.50;
   case 24
      idx = find(jdate - fix(jdate)==0);
      return;
   otherwise
      error('Only H = 6, 12, or 24 is currently supported!');
end

list = mod(jdate,jh);  % modulus
idx  = find(list == 0);


