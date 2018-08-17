function f_rename(x,y)
% - rename variable without memory reallocation
% 
% USAGE: f_rename OldVar NewVar
% 
% This function renames the variable named 'OldVar' to 'NewVar'
% If NewVar already exists, it is overwritten.
% 
% SEE ALSO: f_renameField

% -----References:-----
% This function is a slight modification of D.C. Hanselman's
% <MasteringMatlab@yahoo.com> RENVAR function, with an addition by Tao Lan to
% make sure you don't clear your variable by specifying the same name for both
% old and new names.

% -----Author:-----
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

if nargin~=2 || ~ischar(x) || ~ischar(y)
   error('Two String Arguments Required.')
end
try
   if ~strcmp(x,y) % added by Tao Lan
      evalin('caller',[y '=' x '; clear ' x])
   end
catch
   error('Variable Rename Failed.')
end
