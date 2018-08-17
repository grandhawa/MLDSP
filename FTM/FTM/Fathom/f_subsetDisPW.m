function [sDis,sX,pairList] = f_subsetDisPW(yDis,x)
% - extract subsets of distance matrix based on all pairs of a grouping factor
% 
% USAGE: [sDis,sX,pairList] = f_subsetDis(dis,x);
% 
% Ydis = symmetric distance matrix
% x    = column vector of integers specifying treatment levels
%        or grouping factor
%
% sDis     = cell array of subsets of yDis corresponding to
%            all pairwise treatment levels specified in x
% sX       = cell array of pairwise treatment levels
% pairList = list of pairwise treatement levels
%
% SEE ALSO: f_npManovaPW, f_npManova, f_anosim, f_anosimSub, f_subsetPW

% -----Notes:-----
% This function is used to extract all portions of a square symmetric distance
% matrix based on all pair-wise combinations of treatments levels (or grouping
% factors) specified in the input column vector. It was primarily written as a
% utility function for f_npManovaPW.

% -----Author:-----
% by David L. Jones, Nov-2002
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Mar-2010: pairList is sorted by rows of BOTH columns now; preallocate sDis &
%           sX
% Nov-2011: recode x to avoid errors in f_npManovaPW

% -----Check input:-----
if (size(x,2)>1)
	error('X must be a column vector!');
end

if (size(yDis,1) ~= size(x,1))
	error('yDis & X need same # of rows');
end;

if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix');
end

% Recode grouping factor as consecutive integers:
x = f_recode(x);
% ----------------------

noGrps   = length(unique(x));  % get # of unique treatment levels
pairList = combnk(1:noGrps,2); % get list of all pairwise combinations
pairList = sortrows(pairList,[1 2]); 
noPairs  = size(pairList,1);   % # of pairwise combinations

% Extract portions of yDis corresponding to each pair of treatment levels:
sDis{noPairs} = NaN; % preallocate
sX{noPairs}   = NaN; 
for i = 1:noPairs
	index = find((x == pairList(i,1)) | (x == pairList(i,2)) ); % subset indices of yDis to extract
	
	% Extract subset of distance matrix corresponding to this pair:
	sDis{i} = yDis(index,:);    % rows to keep
	sDis{i} = sDis{i}(:,index); % columns to keep
	
	% Extract subsets of treatment levels:
	sX{i}   = x(index);
end
