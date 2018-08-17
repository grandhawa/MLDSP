function y = f_smooth(x,span,rep)
% - smooth columns of a matrix using a Moving Average filter
%
% USAGE: y = f_smooth(x,span,rep);
%
% x    = matrix 
% span = size of filter (default = 5), should be odd
% rep  = # of times to repeat filter
%
% y   = smoothed data

% -----Notes:-----
% This function is used to smooth columns of a matrix via a Moving Average.
% It utilizes the Matlab function SMOOTH from the Curve Fitting Toolbox.
% Each column of data is smoothed separately but using the same span.
%
% Blocks of data separated by rows of NaN's are smoothed SEPARATELY.
%
% The SMOOTH function differs from similar filters, such as ones that utilize
% the Matlab function FILTER, as it PRESERVES the endpoints of the data set.
%
% Using the repeat option (REP) to run a smaller-sized filter multiple times
% may provide better results than running a larger-sized filter once.

% -----Dependencies:-----
% This function requires the Curve Fitting Toolbox.

% -----Author(s):-----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2003: fixed error with 'remove trailing NaN's'


% -----Check input & set defaults:-----
if (nargin < 2), span = 5; end; % default span size
if (nargin < 3), rep  = 1; end; % don't repeat by default

if (f_isOdd(span)==0)
   error('SPAN must be an ODD number!');
end
% -------------------------------------

[nRows,nCols] = size(x);

% -----Process first & last rows:-----
% remove trailing NaN's:
if (isnan(x(nRows,1))==1)
   x(nRows,:)=[];
   delNaN = 1; % remember we deleted this
else
   delNaN = 0;
end

% make sure first row is NaN's
if (isnan(x(1,1))==0)
   x = [repmat([NaN],1,nCols);x];
   addedNaN = 1; % remember we added this
else
   addedNaN = 0;
end
% ------------------------------------

% get row indices of NaN's:
k = find(isnan(x(:,1))==1);
k = k(:);

% Get indices of blocks of data separated by NaN's:
k2  = k-1;       
k2  = [k2(2:end);size(x,1)];
idx = [k k2]; % indices specifying [start end] of blocks

noBlocks = size(idx,1); % # of blocks

% Extract each block separately into a cell array:
for i=1:noBlocks
   blocks{i} = x(idx(i,1)+1:idx(i,2),:);
end

% -----Smooth each block separately:-----
y = repmat([NaN],1,nCols); % initialize with row of NaN's
for j=1:noBlocks
   subY = subSmooth(blocks{j},span,rep);
   % separate rows blocks with row of NaN's:
   y = [y;repmat([NaN],1,nCols);subY];
end
% ---------------------------------------

% Return to original state:
if (addedNaN==1)
   % strip row of NaN's added to top:
   y = y(2:end,:);
end

if (delNaN==1)
   % replace row of NaN's deleted bottom:
   x(end,:)==repmat([NaN],1,nCols);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTIONS    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = subSmooth(x,span,rep)

nCols = size(x,2);

y = x;

% -----smooth and repeat:-----
for j = 1:rep
   for i=1:nCols
      y(:,i) = smooth(y(:,i),span);
   end     
end
