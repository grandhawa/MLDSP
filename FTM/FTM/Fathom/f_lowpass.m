function Xs = f_lowpass(X,freq)
% - smooth columns of a matrix via a lowpass filter
%
% USAGE: Xs = f_lowpass(X,freq);
%
% X    = input matrix (rows = obs, cols = variables)
% freq = data frequency
%           (e.g., freq = 24 for samples taken every hour)
%           (e.g., freq = 72 for samples every 20 min)
%
% Xs   = smoothed data
% 
% SEE ALSO: f_lowpass40

% -----Notes:-----
% This function uses a Butterworth filter to smooth times-series
% data. It was primarily written to smooth wind and current-meter velocity
% vectors in order to remove events that occur > twice a day, such as
% tides, etc.

% -----References:-----
% after Matlab code by Cynthia Yeung (NOAA-NMFS)

% -----Dependencies:-----
% This function requires the Matlab Signal Processing Toolbox.

% -----Author:-----
% by David L. Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Create Butterworth digital filter:
[b,a] = butter(1,1/(freq/2));

% This is the same as (see f_lowpass40):
% [b,a] = butter(1,2*(1/24)); => this is a 24-hr, 1st order filter

[nr,nc] = size(X);
Xs      = zeros(nr,nc); % preallocate results array

% Apply digital filter:
for i=1:nc
   Xs(:,i) = filtfilt(b,a,X(:,i));
end


