function [sY,sX,pairList] = f_subsetPW(Y,X)
% - extract subsets of column vector based on all pairs of a grouping factor
% 
% USAGE: [sY,sX,pairList] = f_subsetDis(Y,X);
% 
% Y = column vector of response variable
% X = column vector of whole numbers specifying group membership
%
% sY       = cell array of subsets of Y corresponding to
%            all pairwise treatment levels specified in X
% sX       = cell array of pairwise treatment levels
% pairList = list of pairwise treatement levels
%
% SEE ALSO: f_subsetDisPW, f_ancovaPW

% -----Notes:-----
% This function is used to extract all portions of a column vector based on
% all pair-wise combinations of treatments levels (or grouping factors)
% specified in the input column vector. It was primarily written as a 
% utility function for f_ancovaPW.

% -----Author(s):-----
% by David L. Jones, Nov-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (size(Y,2)>1) || (size(X,2)>1)
	error('X and Y must be column vectors!');
end

if (size(Y,1) ~= size(X,1))
	error('Y & X need same # of rows');
end;

% Recode grouping factor as consecutive integers:
X = f_recode(X);
% ----------------------

noGrp    = length(unique(X));  % get # of unique treatment levels
pairList = combnk(1:noGrp,2); % get list of all pairwise combinations
pairList = sortrows(pairList,[1 2]); 
noPairs  = size(pairList,1);   % # of pairwise combinations

% Extract portions of Y corresponding to each pair of treatment levels:
sY{noPairs} = NaN; % preallocate
sX{noPairs} = NaN; 
for i = 1:noPairs
	idx = find((X == pairList(i,1)) | (X == pairList(i,2)) ); % subset indices of Y to extract
	
	% Extract subset of distance matrix corresponding to this pair:
	sY{i} = Y(idx); % rows to keep
		
	% Extract subsets of treatment levels:
	sX{i}   = X(idx);
end
