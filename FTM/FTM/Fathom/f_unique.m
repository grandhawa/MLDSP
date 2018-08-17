function U = f_unique(X)
% - returns unsorted list of unique values
% 
% USAGE: U = f_unique(X)
% 
% X = iput matrix
% U = unsorted list of unique values in X
% 
% SEE ALSO: unique

% -----References:-----
% http://www.mathworks.de/matlabcentral/newsreader/view_thread/236866

% -----Author:-----
% by David L. Jones, Sep-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jan-2011: added support for cell arrays

if iscell(X)
   [nul,I] = unique(X,'first');
else
   [nul,I] = unique(X,'rows','first');
end

U = X(sort(I),:);
