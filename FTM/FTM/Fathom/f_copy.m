function xString = f_copy(x,tab)
% - copies a numerical-array to the clipboard
%   
% USAGE: xString = f_copy(x,tab);
% 
% x   = numerical array to copy to the clipboard as a space (default) or tab
%       separated string. Space-delimited is suitable for pasting into a text
%       editor, while tab-delimited is good for spreadsheets.
% tab = use TAB as the delimiter (default = 0)
% 
% xString = space- or tab-delimited result included for completeness.
% 
% SEE ALSO: f_table

% -----References:-----
% Slightly modified after Matlab function 'num2clip' by Grigor Browning,
% 02-Sept-2005 

% -----Author:-----
% by David L. Jones, Jan-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), tab =  0; end % defalut space-delimit

% convert the numerical array to a string array
% note that num2str pads the output array with space characters to account
% for differing numbers of digits in each index entry
xString = num2str(x); 

% add a carrige return to the end of each row:
xString(:,end+1) = char(10);

% reshape the array to a single line, note that the reshape function reshape is
% column based so to reshape by rows one must use the inverse of the matrix;
% reshape the array to a single line
xString = reshape(xString',1,numel(xString));

% create a copy of arraystring shifted right by one space character:
xStringShift = [' ',xString];

% add a space to the end of arraystring to make it the same length as
% arraystringshift:
xString = [xString,' '];

% now remove the additional space charaters - keeping a single space
% charater after each 'numerical' entry
xString = xString((double(xString)~=32 | double(xStringShift)~=32) &...
   ~(double(xStringShift==10) & double(xString)==32) );

if tab>0
   % convert the space characters to tab characters:
   xString(double(xString)==32) = char(9);
end

% copy the result to the clipboard ready for pasting:
clipboard('copy', xString);
