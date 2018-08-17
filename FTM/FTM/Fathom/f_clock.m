function result = f_clock
% - current date and time as hyphen ('-') delimited character string
% 
% USAGE: result = f_clock
% 
% result = character string in 'yyyy-mo-dd-hh-mm-ss' format
% 
% SEE ALSO: clock, f_renameBatch

% -----Author:-----
% by David L. Jones, Aug-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Get date, time:
T = clock;

% Parse into separate fields:
YY = sprintf('%02.0f',T(1));
MO = sprintf('%02.0f',T(2));
DD = sprintf('%02.0f',T(3));
HH = sprintf('%02.0f',T(4));
MM = sprintf('%02.0f',T(5));
SS = sprintf('%02.0f',T(6));

% Wrap results up into a character string:
result = [YY '-' MO '-' DD '-' HH '-' MM '-' SS];
