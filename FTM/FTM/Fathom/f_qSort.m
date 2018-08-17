function [X,idx] = f_qSort(X,i,j)
% - quick sort elements i:j of a vector, ascending
% 
% USAGE: [sX,idx] = f_qsort(X,i,j);
% 
% X = input vector to sort
% i = initial element of X to sort
% j = final element of X to sort
% 
% sX  = sorted version of vector X
% idx = sorted index to original elements

% -----Notes:-----
% This function was written to replicate the functionality of the R function
% 'R_qsort_I' which is provided as 'qsort.c' in the randomForest package for R.

% -----Author:-----
% by David L. Jones, Jan-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
n = numel(X); % get size of X

if isvector(X)==0, error('X must be a column vector!'); end
if j<i, error('I,J must indicate beginning:ending elements to sort!'); end
if j>n, error('J exceeds the number of elements of X!'); end

% Sort:
idx     = (1:n)';                       % create index to original sort order
sX       = sortrows([X(i:j) idx(i:j)]); % sort
X(i:j)   = sX(:,1);                     % insert sorted portion of X
idx(i:j) = sX(:,2);                     % insert sorted portion of idx
