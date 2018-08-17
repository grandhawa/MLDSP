function idx = f_findCell(x,y)
% - index to rows of Y that match X
%
% USAGE: idx = f_findCell(X,Y)
%
% X = cell array to look for ( = the pattern)
% Y = cell array to look in  ( = array searched)
%
% idx = index to rows of Y that match X

% -----Author:-----
% by David L. Jones, Jun-2007
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% June-2007: now uses intersect (vs. ismember)
% Aug-2010:  revert back to ismember as intersect sorts the result

idx = find(ismember(y,x) == 1);
% [C,idx] = intersect(y,x);
% [C,idx] = intersect(x,y);