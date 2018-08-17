function result = f_numAtLength(num, aveTL, minTL, maxTL, adj, plt)
% - estimate number of fish within length classes for visual survey data
%
% USAGE: result = f_numAtLength(num, aveTL, minTL, maxTL, adj, plt)
%
% num   = number of fish at collection site (column vector)
% aveTL = corresponding average length
% minTL = corresponding minimum length
% maxTL = corresponding maximum length
% adj   = adjust abundances to whole numbers (default=0)
% plt   = plot data for first row            (default=0)
%
% result = output structure with the following fields:
%  .L   = length-classes (up to 5)
%  .A   = corresponding abundances

% -----Notes:-----
% This function estimates the number of fish within up to five length classes
% for visual survey data where only the total number of fish and their
% minimum, average, and maximum lengths are available. This technique
% employs a triangular distribution-based smoothing function to map these
% data as parameters of a frequency distribution having a mode and lower/upper
% bounds defined by the observed mean size and minimum/maximum sizes,
% respectively.

% -----References:-----
% Ported to Matlab from original SAS code 'number-at-length' by Steve Smith
% <steve.smith@rsmas.miami.edu>
% 
% Jones, D. L., J. F. Walter, E. N. Brooks, and J. E. Serafy. 2010.
%   Connectivity through ontogeny: fish population linkages among mangrove
%   and coral reef habitats. Mar. Ecol. Prog. Ser. 401: 245-258.

% -----Author:-----
% by David L. Jones, Jan-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 5), adj = 0; end; % no adjustment to abundance
if (nargin < 6), plt = 0; end; % no plots by default

% Force input as column vectors:
num   = num(:);
aveTL = aveTL(:);
minTL = minTL(:);
maxTL = maxTL(:);

% Check sizes of input:
if numel(unique([numel(num) numel(aveTL) numel(minTL) numel(maxTL)]))>1
   error('NUM, AVETL, MINTL, MAXTL must all be column vectors of the same size!')
end
% -------------------------------------

n        = size(num,1);
result.L = zeros(n,5); % preallocate
result.A = zeros(n,5);

% Process each taxon-transect combination separately:
for  i=1:n
   
   if num(i)==0
      L = [0 0 0 0 0];
      A = [0 0 0 0 0];
      
   elseif (size(unique([aveTL(i) minTL(i) maxTL(i)]),2)==1) % sizes duplicated
      L = [aveTL(i) 0 0 0 0];
      A = [num(i) 0 0 0 0];
      
   elseif num(i)==1
      L = [aveTL(i) 0 0 0 0];
      A = [1 0 0 0 0];
      
   elseif num(i)==2
      L = [minTL(i) maxTL(i) 0 0 0];
      A = [1 1 0 0 0];
      
   elseif num(i)==3
      L = [minTL(i) aveTL(i) maxTL(i) 0 0];
      A = [1 1 1 0 0];
      
   elseif num(i)==4
      L = [minTL(i) aveTL(i) maxTL(i) 0 0];
      A = [1 2 1 0 0];
      
   elseif num(i)>4
      L = [minTL(i) [minTL(i)+aveTL(i)]/2 aveTL(i) [aveTL(i)+maxTL(i)]/2 maxTL(i)];
      
      % if num(i)<=99
      A(1) = 1;
      A(3) = num(i)/2;
      A(5) = 1;
      A(2) = (num(i)-[A(1)+A(3)+A(5)])/2;
      A(4) = A(2);
      % else % num(i)>99
      % A(1) = 0.01*num(i);
      % A(3) = num(i)/2;
      % A(5) = 0.01*num(i);
      % A(2) = (num(i)-[A(1)+A(3)+A(5)])/2;
      % A(4) = A(2);
      % end
      
      % Make adjustments so abundances are only whole numbers:
      if adj>0
         frac = A(2) - floor(A(2)); % extraction fractional component
         if frac>0
            A(2) = ceil(A(2));        % round up
            A(4) = A(2);              % round up
            A(3) = A(3)-([1-frac]*2); % compensate for the rounding ups
         end
      end
      
   else
      error(['There is a problem with the abundance of ROW ' num2str(i) '!']);
   end
   
   result.L(i,:) = L;
   result.A(i,:) = A;
   
end

if plt>0 % Optionally plot data from row 1:
   figure;
   plot(result.L(1,:)', result.A(1,:)', 'b.-');
   xlabel('Length');
   ylabel('Abundance');
   title('Abundance-at-Length');
end







