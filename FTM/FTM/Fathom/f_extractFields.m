function [result,labels] = f_extractFields(X,idx)
% - extract structure fields & combine into a single matrix
%
% USAGE: [result,labels] = f_extractFields(X,idx);
%
% X      = input structure
% idx    = vector of indices of fields to extract
%          (default = extract all)
%
% result = matrix of data from extracted fields
% labels = cell array of field names
%
% SEE ALSO: f_struct2flat

% -----Notes:-----
% This function was written primarily to extract and combine (column-wise)
% the data in a structure returned from QUERYBUILDER, part of the Matlab 
% Database Toolbox.
%
% Extracted fields are combined column-wise in a matrix, so they must all
% be of the SAME type and size.

% -----Example:-----
% Suppose myData is a structure with 10 fields and you want to EXTRACT and
% COMBINE the last 4 fields into a single matrix:
%
% >> myData
% myDdata = 
%          Site: {452x1 cell}
%          Fish: {452x1 cell}
%     Collected: {452x1 cell}
%         SL_mm: [452x1 double]
%             E: [452x1 double]
%             N: [452x1 double]
%            Li: [452x1 double]
%            Na: [452x1 double]
%            Mg: [452x1 double]
%             P: [452x1 double]
%
% >> [result,labels] = f_extractFields(myData,[7:10]);
%
% Extract and combine the first 3 fields:
%
% >> [result,labels] = f_extractFields(myData,[1:3]);

% -----Author:-----
% by David L. Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2003: updated documentation, set default value for idx

% -----Check input and set defaults:-----
if (isstruct(X)~=1)
   error('X must be a structure!');
end

% Get all variable labels:
all_labels = fieldnames(X);

if (nargin < 2)
   idx = [1:size(all_labels,1)];
end
% ----------------------

idx = idx(:); % force column vector

% Check input:
if [max(idx) > size(all_labels,1)]
   error('IDX exceeds the # of fields in X!');
end

% Keep only required variable labels:
labels = all_labels(idx);

% Build command string to evaluate:
nr        = size(labels,1);
cmdString = [];
for i=1:nr
   cmdString = [cmdString ' X.' labels{i}];
end
cmdString    = ['[' cmdString '];'];
cmdString(2) = []; % remove extra space

% Evaluate command string:
result = eval(cmdString);

