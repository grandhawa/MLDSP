function  Y = f_cellstr2num(C,S,N)
% - convert a cell array of strings to a single categorical variable
%
% USAGE: Y = f_cellstr2num(C,S,N)
%
% C = cell array of strings
% S = unique strings that occur in C
% N = corresponding numbers that values of S should be converted to in C
%
% Y = column vector of integers specifying factor levels or group membership
%
% SEE ALSO: f_cellstr2cat

% -----Author:-----
% by David L. Jones, Jun-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Notes:-----
% The purpose of this code is to create a categorical variable from a cell
% array of strings representing some type of text labeling system). This
% can then be used as a grouping variable in MANOVA, f_dummy, f_grpMean, etc.
%
% Ths function is similar to f_cellstr2cat, but allows you to explicitly
% determine the resulting number that corresponds to each unique string

% -----Check input:-----
if ~iscell(C), error('C must be a cell array of strings!'); end
if ~iscell(S), error('S must be a cell array of strings!'); end

if numel(S) ~= numel(N)
   error('S and N must have the same number of elements!');
end
% ----------------------

C  = C(:);      % force column vector
nr = size(C,1); % get # rows
nG = numel(S);  % # groups
Y  = nan(nr,1); % preallocate

for i=1:nG
   idx    = ismember(C,S(i)); % index to rows of C belonging to this group in S
   Y(idx) = N(i);             % convert to corresponding number
end

% Make sure all strings in C were replaced:
if sum(isnan(Y))>0
   error('There were string in C not specified by S!')
end
