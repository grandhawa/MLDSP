function [Y,S,res] = f_rosner(X,ave,s,p,crit)
% - Rosner's test for spike elimination (many outlier removal)
%
% USAGE: [Y,S,res] = f_rosner(X,ave,s,p,crit);
%
% X    = input data matrix (rows = observations, cols = variables)
% ave  = replace spikes with smoothed data (= 1) or NaN (= 0)         (default = 1)  
% s    = strength of smoothing filter on a relative scale from 0 to 1 (default = 0.5)
% p    = # of outliers to evaluate, as a proportion of sample size    (default = 0.1)
% crit = detection criterion lambda                                   (default = 2)
%
% Y   = output data with spikes (outliers) replaced with moving average or NaN
% S   = logical array where 1 indicates original datum was an outlier
% res = residuals
%
% SEE ALSO: f_grubbs, f_spike

% -----References:-----
% Olsson, G., and B. Newell. 1999. Wastewater treatment systems: modelling,
% diagnosis and control. IWA Publishing, London. pp.140-141
% 
% Matlab code after outlier.m by Bob Newell, February 1996

% -----Author:-----
% by David L. Jones, Oct-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2011: replaces spikes with 'smoothed' value vs. overall mean, now
%           outputs residuals

% -----Check input & set defaults:-----
if (nargin < 2), ave  =  1;   end % default replace outliers smoothed data
if (nargin < 3), s    =  0.5; end % default strength = 0.5
if (nargin < 4), p    =  0.1; end % default # of outliers to evaluate 
if (nargin < 5), crit =  2;   end % default crit of 2

% Check proportions:
if (p<=0) || (p>1)
   error('P must be between 0 & 1!');
end
% -------------------------------------

[r,c] = size(X);            % get size of input
S     = zeros(r,c);         % initialize logical array

% Filter input data:
% Xs  = f_filterMA(X,span,1); % smoothed via a moving average
Xs  = f_butter(X,8,s);        % light lowpass Butterworth
res = X - Xs;                 % residuals

k   = ceil(r*p);              % # of outliers to evaluate

% Find/Replace outliers for each column separately:
for i=1:c
   R = zeros(k+1,1);
   
   % Sort deviations from the mean:
   xBar    = mean(res(:,i));
   [xs,is] = sort(abs(res(:,i) - xBar));
   
   % Calculate statistics for up to k outliers:
   for j = 0:k,
      xx     = xs(1:r-j);
      R(j+1) = abs( xx(r-j) - mean(xx) ) / std(xx);
   end
   
   % Statistical test to find outliers:
   idx = [];
   for j = 1:k,
      if R(j)/R(k+1) > crit,
         idx = [idx is(r-j+1)];
      end
   end
   
   % Replace outliers:
   if (~isempty(idx))
      if (ave>0)
         X(idx,i) = Xs(idx,i); % replace outlier with moving average
      else
         X(idx,i) = NaN;       % replace outlier with missing value
      end
      S(idx,i) = 1;            % mark values that were replaced
   end
end

% Rename for output:
Y = X;
