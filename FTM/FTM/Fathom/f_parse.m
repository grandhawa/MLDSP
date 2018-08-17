function result = f_parse(x,labels)
% - parse matrix to structure
% 
% USAGE: result = f_parse(x,labels);
% 
%  x     = input matrix (rows = obs, col = variables)
% labels = cell array of variable labels
%           e.g., labels = {'sal' 'tmp' 'elev'};

% -----Notes:-----
% labels cannot begin with a number

% -----Author:-----
% by David L. Jones, Oct-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----

% if labels are not cell arrays, try forcing them:
if ~iscell(labels), error('LABELS must be a cell array!'); end;

labels = labels(:); % force as row vector
if (size(x,2) == size(labels,1))<1
   error('The # of rows in X & LABELS must be equal!');
end;
% -----------------------------------------------

nCol = size(x,2); % get # columns

for i=1:nCol
   cmdStr = sprintf('result.%s = x(:,%d);',labels{i},i); % command string
   eval(cmdStr);                                        % evaluate command
end