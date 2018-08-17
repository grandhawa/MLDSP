function [Y,S] = f_grubbs(X,ave,alpha,span)
% - Grubbs test for spike elimination (outlier removal)
%
% USAGE: [Y,S] = f_grubbs(X,ave,alpha,span);
%
% X     = input data matrix (rows = observations, cols = variables)
% ave   = replace spikes with a moving average (= 1) or NaN (= 0)  (default = 1)
% alpha = significance level                                       (default = 0.05)
% span  = size of moving average window used when ave=1            (default = 5)
%
% Y = output data with spikes (outliers) replaced by a moving average or NaN
% S = logical array where 1 indicates original datum was an outlier
%
% SEE ALSO: f_spike, f_rosner

% -----References:-----
% modified after Brett Shoelson's 'deleteoutliers.m'
%
% Copyright (c) 2009, Brett Shoelson
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% -----Author:-----
% modifications by David L. Jones, Aug-2010
%
% This file is part of the FATHOM Toolbox for Matlab.

% Apr-2011: outliers now replaced by moving average (vs. overall average);
%           reorganized code

% -----Check input & set defaults:-----
if (nargin < 2), ave   =  1;    end % default replace outliers wih moving average
if (nargin < 3), alpha =  0.05; end % default alpha of 0.05
if (nargin < 4), span  =  5;    end % default span of 5
% -------------------------------------

[r,c] = size(X);
S     = zeros(r,c); % initialize logical array

% Critical value for Grubbs Test:
tcrit   = tinv(alpha/(2*r),r-2);
critVal = (r-1)/sqrt(r)*(sqrt(tcrit^2/(r-2+tcrit^2)));

if (ave>0)
   Xs = f_filterMA(X,span,1); % moving average
end

% Find/Replace outliers for each column separately:
for i=1:c
   outlier = 1;
   while outlier
      Xc      = X(:,i);                                        % extract 1 column of X
      aveVar  = nanmean(Xc);
      idx     = find( abs(Xc-aveVar) == max(abs(Xc-aveVar)) ); % find largest spike
      maxVar  = Xc(idx(1));                                    % in case of ties, choose one
      sdVar   = std(Xc);
      tn      = abs((maxVar-aveVar)/sdVar);
      
      outlier = (tn > critVal);
      
      if outlier
         X(idx,i) = NaN; % mark as missing
         S(idx,i) = 1;   % update index to outliers
      end
   end
end

Y = X; % make copy for output

% Optionally replace outlier with moving average:
if (ave>0)
   for i=1:c
      idx      = find(isnan(X(:,i)));
      Y(idx,i) = Xs(idx,i);
   end
end
