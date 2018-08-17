function idx = f_sub(jdate)
% - subsample data every 6 hrs
%
% USAGE: idx = f_sub(jdate);
%
% x   = vector of julian dates
% idx = index to elements that occur every 6 hrs.
%
% SEE ALSO: f_julian

% -----Author:-----
% David L. Jones, Oct-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

idx = find(mod(jdate,0.25) == 0);

  