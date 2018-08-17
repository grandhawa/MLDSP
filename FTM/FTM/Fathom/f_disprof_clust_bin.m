function result = f_disprof_clust_bin(A,B)
% - compare 2 DISPROF binary connectivity matrices
%
% USAGE: result = f_disprof_clust_bin(A,B);
%
% A,B = binary connectivity matrices (BIN) created by f_disprof_clust
%
% result = structure of results with the following fields:
%   .n   = # of objects in A that differ from B
%   .idx = index to rows of A that differ from B
%   .cng = measure of congruence between A & B (ranges from 0-1)
%
% SEE ALSO: f_disprof_clust

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
% Check sizes of input:
if ~isequal(size(A),size(B))
   error('A & B are not of compatible sizes!')
end

% Make sure input is a binary connectivity matrix:
if (f_issymdis(A) == 0)
   error('A & B must be square symmetric distance matrices!');
end

% Check that inputs are binary:
if any(((A==0)+(A==1)<1) | ((B==0)+(B==1)<1))
   error('A & B must consist of only 1''s or 0''s!')
end

% Check if inputs are identical:
if isequal(A,B)
   result.n   = 0;
   result.idx = [];
   result.cng = 1;
   return
end
% -------------------------------------

% Measure the congruence between A & B:
X     = double(A~=B);                    % find elements that do not match
ai    = logical(tril(ones(size(X)),-1)); % get indices for lower diagonal
X(ai) = NaN;                             % replace lower tridiagonal with NaN's
idx   = find(nansum(X,2));               % get index to rows that differ
n     = numel(idx);                      % get of rows that differ
cng   = 1 - n/size(A,1);                 % get proportion of agreement

% Wrap results up into a structure:
result.n   = n;
result.idx = idx;
result.cng = cng;
