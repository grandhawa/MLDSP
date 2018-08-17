function Xs = f_butter(X,n,s)
% - smooth columns of time series data via Butterworth lowpass filter
%
% USAGE: Xs = f_butter(X,n,s);
%
% X = input matrix (rows = obs, cols = variables)
% n = Nth order of lowpass filter
% s = strength of smoothing filter on a relative scale from 0 to 1 (default = 1)
%
% Xs = smoothed data
% 
% SEE ALSO: f_lowpass, f_lowpass40, f_filterMA

% -----Notes:-----
% This is a generic smoothing function that filters time series data with a
% lowpass Butterworth filter equivalent to a 40-h lowpass filter when data are
% sampled hourly (and there are 24 hrs/day). However, the frequency of the
% sampling interval of the input data are not required in this function, only
% that we use a similar ratio (i.e., 1/40).
% 
% The input variable 's' allow setting the relative strength of the filter, with
% values of 1 providing a full smoothing effect equivalent to a 40-h lowpass
% filter.

% -----References:-----
% - email/Matlab code by Bill Johns<bjohns@rsmas.miami.edu>
% - filt.m in Olsson, G. and Newell. 1999. Wastewater Treatment Systems:
%   Modelling, Diagnosis and Control

% -----Dependencies:-----
% This function requires the Matlab Signal Processing Toolbox.

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 3), s =  1; end % default full strength

if (s<=0) || (s>1)
   error('S must be greater than 0 but not 1!');
end

[nr,nc] = size(X);
Xs      = zeros(nr,nc); % preallocate results array

% Create Butterworth digital filter equivalent to a 40-h lowpass filter when
% data are sampled hourly:
[b,a] = butter(n,2*(1/(40*s))); 

% Apply digital filter:
for i=1:nc
   Xs(:,i) = filtfilt(b,a,X(:,i));
end
