function [Xs,idxS,idxU] = f_sort(X)
% - sort elements of X (ascending) and create sorting/unsorting indices
% 
% USAGE: [Xs,idxS,idxU] = f_sort(X);
% 
% X    = input vector or matrix (columns sorted separately)
% 
% Xs   = sorted version of X
% idxS = index to repeat the sort; i.e.,  X(idxS) = Xs
% idxU = index to undo the sort;   i.e., Xs(idxU) = X
% 
% SEE ALSO: sort, sortrows

% -----References:-----
% http://blogs.mathworks.com/loren/2007/08/21/reversal-of-a-sort/

% -----Author:-----
% by David L. Jones, Mar-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Sort X and return sorting index:
[Xs,idxS] = sort(X);

% Create unsorting index:
idxU(idxS) = 1:size(X,1);
