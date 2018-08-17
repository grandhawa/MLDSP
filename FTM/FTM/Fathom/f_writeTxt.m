function f_writeTxt(x,fname)
% - write character array to a text file
%
% USAGE: f_writeTxt(x,'fname');
%
% x     = char array
% fname = name of file to write to

% -----Author:-----
% by David L. Jones, Aug-2007
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (ischar(x) == 0)
   error('X must be a CHAR array!');
end

fid = fopen(fname, 'w');
fwrite(fid, x, 'char');
fclose(fid);

