function Xs = f_lowpass40(X,freq,M)
% - smooth columns of a matrix via a 6th order, 40-hr lowpass filter
%
% USAGE: Xs = f_lowpass40(X,freq,M);
%
% X    = input matrix (rows = obs, cols = variables)
% freq = sampling interval
%           (e.g., freq =   1 for samples taken every hour)
%           (e.g., freq = 0.5 for samples taken every half hour)
% M    = value to replace missing (= NaN) data with (default = 0)
%
% Xs = smoothed data
% 
% SEE ALSO: f_lowpass

% -----Notes:-----
% This function uses a Butterworth filter to smooth times-series data. It was
% primarily written to smooth wind and current-meter velocity vectors in order
% to remove events that occur > twice a day, such as tides, etc.
% 
% Dr. Johns usually uses a 4th or 6th order Butterworth filter passed forward
% and backward over the data using 'filtfilt' to be sure to eliminate phase
% shifts.
% 
% When NaN's are present in X, the filter returns only NaN's, so you must
% replace them temporarily something, run the filter, then return the NaN's
% appropriately.

% -----References:-----
% after email/Matlab code by Bill Johns<bjohns@rsmas.miami.edu>

% -----Dependencies:-----
% This function requires the Matlab Signal Processing Toolbox.

% -----Author:-----
% by David L. Jones, Dec-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check inputs and set defaults:-----
if (nargin < 3), M = 0; end % replace NaN's with 0 by default
if size(X,2)>1, error('X must be a column vector!'); end

% Replace Nan's if they are present:
if sum(isnan(X))>0,
   NaNflag   = 1;                 % set a flag that there were missing data
   idxNaN    = find(isnan(X)==1); % get indices to missing data
   X(idxNaN) = M;                 % replace with M
else
   NaNflag = 0;
end

order  = 4;  % 6th order
cutoff = 40; % 40 hour
% ----------------------------------------

% Create Butterworth digital filter:
[b,a] = butter(order,2*freq/cutoff);

[nr,nc] = size(X);
Xs      = zeros(nr,nc); % preallocate results array

% Apply digital filter:
for i=1:nc
   Xs(:,i) = filtfilt(b,a,X(:,i));
end

% Return NaN's to data as appropriate:
if (NaNflag>0)
   Xs(idxNaN) = NaN;
end
   
   
