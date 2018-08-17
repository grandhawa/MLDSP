function X = f_recode(Y)
% - recode elements of vector as consecutive integers
% 
% USAGE X = f_recode(Y);
% 
% Y = vector or matrix of input variables (column-wise) 
% X = vector or matrix recoded as consecutive integers
%
% SEE ALSO: f_designMatix, dummvar

% -----Notes:-----
% This function is primarily a fix for f_designMatrix which uses Matlab's
% DUMMYVAR function. It is used, for example, to convert treatment levels [2 2 2
% 5 5 5] to [1 1 1 2 2 2] so proper ANOVA design matrices can be constructed.

% -----Author:-----
% by David L. Jones, Nov-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2008: removed 'sort' from 'unique' command as it's alread sorted; use
%           logical indexing vs. find; input/output are treated as separate
%           variables; added support for matrices

[r,c] = size(Y);
X = repmat(NaN,r,c); % preallocate

for i = 1:c
   X(:,i) = sub_recode(Y(:,i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = sub_recode(Y)
Y  = Y(:);      % make sure it's a column vector
G  = unique(Y); % get sorted list of unique values
nr = size(G,1);   

X = repmat(NaN,nr,1); % preallocate
for i = 1:nr
	X((Y == G(i))) = i;
end
