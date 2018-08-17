function idx = f_firstOccur(vector)
% - returns indices of the first occurrence of unique elements of input vector
%
% USAGE: idx = f_firstOccur(vector)

% by David L. Jones, July-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% June-2006: added the sort to keep original order of vector since 'unique'
% function returns a SORTED list of unique elements

theElements = unique(vector); % get unique elements of input vector
theElements = theElements(:); % make sure it's a column vector
noRows = size(theElements,1); 

for i = 1:noRows
   idx(i) = min(find(theElements(i) == vector));
end

idx = sort(idx(:)); % sort and make sure it's a column vector