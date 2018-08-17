function y = f_wtMean(wt,x, method)
% - weighted mean
% 
% USAGE: y = f_wtMean(wt,x,method)
% 
% wt     = column vector of weighting factors (non-negative)
% x      = input data (column vector or matrix)
% method = 1: wts are weights (default)
%          2: wts are integers specifying a 'repeat' or 'effort' factor
% 
% y  = weighted mean (scalar or row vector)

% -----References:-----
% http://en.wikipedia.org/wiki/Weighted_mean

% -----Notes:-----
% Valus with a high weight contribute more to the weighted mean than do elements
% with a low weight.

% -----Author:-----
% by David L. Jones, Feb-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), method   = 1; end; % default method

if size(wt,1) ~= size(x,1)
   error('WT and X must have the SAME # of rows!');
end

if min(wt)<0
   error('Negative weighting factors (WT) are not allowed!')
end
% -------------------------------------

if method == 1

   y = sum(x.*repmat(wt,1,size(x,2)))/sum(wt);

elseif method == 2

   n = size(wt,1); % get # rows
   z = [];         % preallocate

   for i=1:n
      z = [z; repmat(x(i,:),1)];
   end

   y = mean(z);

else
   error('Method must be 1 or 2!');
end
