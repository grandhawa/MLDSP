function cmd = f_renameField(S,x,y)
% - rename a structure field
%
% USAGE: f_renameField(S,'x','y')
%
% S = input structure
% x = character string specifying name of an existing field of S
% y = character string specifying a new field name
%
% SEE ALSO: f_rename

% -----References:-----
% This function is modified after 'f_rename' and ideas from:
% http://blogs.mathworks.com/loren/2010/05/13/rename-a-field-in-a-structure-array/

% -----Author:-----
% by David L. Jones, Aug-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
% Make sure S is a structure:
if (~isstruct(S)==1)
   error('S must be a valid structure!');
end

% Make sure S.x is a field:
if (~isfield(S,x)==1)
   error('X must be a field of structure S!');
end

% Make sure S.y doesn't already exist:
if (isfield(S,y)==1) && (~strcmp(x,y))
   error('S.Y already exists, rename S.X something else!');
end

if nargin < 3 || ~ischar(x) || ~ischar(y)
   error('Two String Arguments Required.')
end
% ---------------------------------------

% Replace S with the name of the input structure:
S = inputname(1);

% Rename structure field:
try
   if (~strcmp(x,y)) % skip if both specify the same name
      % Generate command string:
      cmd = [S '.' y ' = ' S '.' x '; ' S ' = rmfield(' S ',''' x ''');'];
      
      % Evaluate command string in caller's workspace:
      evalin('caller',cmd);
   end
catch
   error('Failed to rename the structure''s field!')
end
