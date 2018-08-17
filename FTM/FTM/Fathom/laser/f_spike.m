function [Y,S] = f_spike(X,ave,n)
% - spike elimination (outlier removal)
% 
% USAGE: [Y,S] = f_spike(X,ave,n);
% 
% X   = input data matrix (rows = observations, cols = variables)
% ave = replace spikes with average values (= 1, default) or NaN (= 0)
% n   = # of standard deviations defining threshold for removal (default = 3)
% 
% Y   = output data with spikes (outliers) replaced by the mean
% S   = logical array where 1 indicates original datum was an outlier
% 
% SEE ALSO: f_grubbs, f_rosner

% -----Author:-----
% by David L. Jones<djones@marine.usf.edu>, Oct-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), n  =  3; end % default threshold to filter set at 3 SD's
% -------------------------------------

[r,c] = size(X);
S     = zeros(r,c); % initialize logical array

stdVar = repmat((std(X)),r,1);             % standard deviation
aveVar = repmat((mean(X)),r,1);            % mean
idx    = find(abs(X - aveVar) > stdVar*n); % values more than n SD's from mean

if (ave>0)
   X(idx) = aveVar(idx);                   % replace outliers with column means
else
   X(idx) = NaN;                           % replace outlier with missing value
end

S(idx) = 1;                                % mark values that were replaced
Y      = X;                                % rename for output
