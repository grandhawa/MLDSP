function f_export(filename, m, dlm)
% - write ASCII delimited file.
%   
% Usage: f_export('FILENAME',M,DLM)
% 
% writes matrix M into FILENAME using the character DLM as the delimiter
% (defaults to space-delimited); Specify '\t' to produce tab-delimited
% files, ',' for comma-delimited.
%
% SEE ALSO: f_import

% modified after dlmwrite.m to export values = 0.

% by David L. Jones, April-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


% test for proper filename:
if ~isstr(filename),
   error('FILENAME must be a string.');
end;

if nargin < 2, error('Requires at least 2 input arguments.'); end

NEWLINE = sprintf('\n');

% delimiter defaults to Space (DLJ)
if nargin < 3, dlm = ' '; end
dlm = sprintf(dlm); % Handles special characters.


% open the file
if strncmp(computer,'MAC',3)
   fid = fopen(filename ,'wt');
else
   fid = fopen(filename ,'wb');
end

if fid == (-1), error(['Could not open file ' filename]); end

% simplify function by not using offsets (DLJ)
r = 0;
c = 0;

% dimensions size of matrix
[br,bc] = size(m);

% start with offsetting row of matrix
for i = 1:r
   for j = 1:bc+c-1
      fwrite(fid, dlm, 'uchar');    
   end
   fwrite(fid, NEWLINE, 'char');
end

% start dumping the array, for now number format float
for i = 1:br
   
   % start with offsetting col of matrix
   for j = 1:c
      fwrite(fid, dlm, 'uchar');    
   end
   
   for j = 1:bc
      %% if (m(i,j) ~= 0) %% <---this line removed to allow export of 0's (DJ)
      str = num2str(m(i,j));
      fwrite(fid, str, 'uchar');    
      if(j < bc)
         fwrite(fid, dlm, 'uchar');    
      end
   end
   fwrite(fid, NEWLINE, 'char'); % this may \r\n for DOS 
end

% close files
fclose(fid);