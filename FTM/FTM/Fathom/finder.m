function finder()
% - open present working directory in a MacOS X's 'FINDER' window
% 
% USAGE: 
% >> finder
% 
% SEE ALSO: saver

% -----Author:-----
% by David L. Jones, May-2007
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2007: use 'dot' notation and eliminate calling EVAL
% Sep-2009: wrap directory in single quote to handle folder names with spaces
% Aug-2012: converted from a script to a function so 'MakeContents.m' will
%           include it in 'Contents.m'

if (ismac==0)
   error('This function is designed for MAC OS X!')
end

% eval(['!open ' pwd]);
eval(['!open ''' pwd ''''])

% !open .;
