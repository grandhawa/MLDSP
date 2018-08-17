function y = f_filterMA(x,n,rep,ave)
% - smooth columns of a matrix using a Moving Average (or Median)
%
% USAGE: y = f_filterMA(x,n,rep,ave);
%
% x    = input matrix (rows = observations, cols = variables)
% n    = span (= # points on either side of current point)         (default = 5)
% rep  = # of times to filter data                                 (default = 1)
% ave  = use average instead of median                             (default = 1)
%
% y    = smoothed data
%
% SEE ALSO: f_filterSinclair, f_butter, f_lowpass, f_lowpass40

% -----Notes:-----
% This function is used to smooth a matrix via a Moving Average. It
% utilizes O. Liungman's SMOOTH function. Each column of data is smoothed
% separately but using the same sized filter.

% -----Author(s):-----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
%
% with help from news://comp.soft-sys.matlab

% Apr-2011: edited documentation; default span = 5
% Jun-2011: updated documentation; added support for median
% Oct-2011: updated documentation; replaced a call to mean when is
%           should've been median 

% -----Check input & set defaults:-----
if (nargin < 2), n   = 5; end % default span = 5
if (nargin < 3), rep = 1; end % don't repeat by default
if (nargin < 4), ave = 1; end % default use average instead of median

% Check ave:
if ~isequal(ave,1) && ~isequal(ave,0)
   error('AVE must be 0 or 1!')
end
% -------------------------------------

nCols = size(x,2);
y     = zeros(size(x)); % preallocate

% Smooth each column separately:
for i=1:nCols
   y(:,i) = subSmooth(x(:,i),n,ave);
end

% -----Repeat filter:-----
if (rep>1)
   for j = 2:rep
      y = f_filterMA(y,n,1,ave);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%             SUBFUNCTIONS                           %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yout = subSmooth(yin,N,ave)
% SMOOTH.M: Smooths vector data.
%			YOUT=SMOOTH(YIN,N) smooths the data in YIN using a running mean
%			over 2*N+1 successive point, N points on each side of the
%			current point. At the ends of the series skewed or one-sided
%			means are used.

%			Olof Liungman, 1997
%			Dept. of Oceanography, Earth Sciences Centre
%			Göteborg University, Sweden
%			E-mail: olof.liungman@oce.gu.se

if (nargin < 2), error('Not enough input arguments!'), end
[rows,cols]       = size(yin);
if min(rows,cols) ~=1, error('Y data must be a vector!'), end
if length(N)      ~=1, error('N must be a scalar!'), end

yin         = (yin(:))';
l           = length(yin);
yout        = zeros(1,l);
temp        = zeros(2*N+1,l-2*N);
temp(N+1,:) = yin(N+1:l-N);

for i = 1:N
   if (ave>0)
      yout(i)     = mean(yin(1:i+N));
      yout(l-i+1) = mean(yin(l-i-N:l));
   else
      yout(i)     = median(yin(1:i+N));
      yout(l-i+1) = median(yin(l-i-N:l));
   end
   
   temp(i,:)     = yin(i:l-2*N+i-1);
   temp(N+i+1,:) = yin(N+i+1:l-N+i);
end

if (ave>0)
   yout(N+1:l-N) = mean(temp);
else
   yout(N+1:l-N) = median(temp);
end

if size(yout)~=[rows,cols]
   yout = yout';
end
