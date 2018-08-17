function num = f_month2num(mon)
% - convert month string to number
%
% USAGE: result = f_month2num(mon)
%
% mon = string indicating month (e.g., 'Jun')
% num = number indicating month

% -----Author:-----
% by David L. Jones <djones@marine.usf.edu>,
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

mon = lower(mon); % convert to lower case
mon = mon(1:3);   % just take first 3 characters

switch mon
   case 'jan'
      num = 1;
   case 'feb'
      num = 2;
   case 'mar'
      num = 3;
   case 'apr'
      num = 4;
   case 'may'
      num = 5;
   case 'jun'
      num = 6;
   case 'jul'
      num = 7;
   case 'aug'
      num = 8;
   case 'sep'
      num = 9;
   case 'oct'
      num = 10;
   case 'nov'
      num = 11;
   case 'dec'
      num = 12;
   otherwise
      error('MON doesn''t specifiy a recognizable month!')
end
