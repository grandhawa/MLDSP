function saver
% - save workspace to file tagged in 'fname' variable to current directory
%
% USAGE:
% >> saver
% 
% SEE ALSO: finder

% -----NOTES:-----
% This scripts requires you to have a filename defined in a string variable
% nammed 'fname', for example:
% 
% > fname = 'myfile';

% -----Author:-----
% by David L. Jones, June-2007
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: check that fname exists
% Aug-2012: converted from a script to a function so 'MakeContents.m' will
%           include it in 'Contents.m' (i.e., added calls to 'evalin')

% -----Check input:-----
if ( evalin('caller',['exist(''fname'')']) == 0 );
   error('The variable ''fname'' specifying the filename is not present!')
else
   if ( evalin('caller',['ischar(''fname'')'])== 0);
      error('fname sould be a character array!');
   end
end
% ----------------------

% Get filename from the workspace: 
fname = evalin('caller','fname');

% Evaluate in the workspace:
evalin('caller',['save ' fname])

% Display where data were saved:
fprintf(['\nFile: ' fname ' has been saved... \n\n']);
