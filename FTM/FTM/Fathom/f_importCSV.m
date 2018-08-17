function result = f_importCSV(fname,labels)
% - import numeric data from a comma-separated-values (*.csv) file
%
% USAGE: result = f_importCSV('fname',labels);
%
% fname   = name of *.csv file to import
% labels  = 1: data file only includes column labels               (default = 1)
%           2: data file includes column and row labels
% 
% result = structure of imported data with the following fields:
%  .dat = numeric data
%  .txt = cell array of strings specifying column labels
%  .row = cell array of strings specifying row labels (the name of this
%         field is from the csv file)
% 
% SEE ALSO: f_importXL

% -----Notes:-----
% The purpose of this function is to import numeric data with column
% (variable) labels into the Matlab workspace as a single, structure
% variable. The column labels are imported as a cell array of strings and
% there is optional support for row labels as well.
% 
% Example *.csv data (with row and col labels) should look like this:
% 
% sites,sp_A,sp_B,sp_C,sp_D,sp_E
% stn_1,11,22,33,44,55
% stn_2,0.2,0.6,0.7,0.1,0.0
% stn_3,122,66,902,7,66
% stn_4,0.9,33,0.04,5,99
% 
% Note: data in an Excel *.xls file can be saved as a comma-separated-values
% (*.csv) file before being imported with this function.

% -----References:-----
% http://blogs.mathworks.com/loren/2010/05/13/rename-a-field-in-a-structure
% -array/

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


% -----Set defaults & check input:-----
if (nargin < 2), labels = 1; end % default only col labels provided

% Check that file 'fname' is in your path:
if (exist(fname,'file')~=2)
   error(['Can''t find file ' fname '!'])
end

% Check labels:
if (labels~=1) && (labels~=2)
   error('LABELS must be 1 or 2!');
end
% -------------------------------------

% Read data:
temp = importdata(fname,',');

% Load numeric data:
result.dat = temp.data;

% Parse row/col labels:
if (labels==1)
   result.txt = temp.textdata(1,1:end);
else
   result.txt = temp.textdata(1,2:end);
   result.row = temp.textdata(2:end,1);
end

% Check for compatibility:
if size(result.dat,2) ~= numel(result.txt)
   fprintf('\nSize of column labels don''t match numeric data! \n\n')
   error('Make sure your value for LABELS is correct!')
end
if (labels==2) && size(result.dat,1) ~= numel(result.row)
   fprintf('\nSize of row labels don''t match numeric data! \n\n')
   error('Make sure your value for LABELS is correct!')
end

% Rename 'row' field:
if (labels==2)
   newField            = temp.textdata{1};
   [result.(newField)] = result.row;
   result              = rmfield(result,'row');
end
