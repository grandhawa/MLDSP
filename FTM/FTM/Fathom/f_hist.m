function h = f_hist(Y)
% - create a normalized histogram plot of a column vector
%
% USAGE: h = f_hist(Y);
%
% Y = column vector of sample observations
% h = graphics handle
%
% SEE ALSO: hist

% -----Notes:-----
% Portions of the plot routine are modified after Mike Sheppard's 'allfitdist'
% function and my 'f_johnson_fit' function.

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
% Check size of input:
if (size(Y,2)>1)
   error('Y must be a column vector!');
end

% Check for missing values:
if any(isnan(Y))
   error('Y contains NaN''s!');
end
% -------------------------------------

colorVar = [1 1 1]*0.65; % define histogram color

figure; set(gcf,'color','w');
hold on;
box  on;
grid on;

% Plot histogram:
nbins = max(min(length(Y)./10,100),50); % get # bins
yi    = linspace(min(Y),max(Y),nbins)'; % range
dy    = mean(diff(yi));                 % step size
yfi   = histc(Y,yi-dy);                 % get frequencies
yfi   = yfi./sum(yfi)./dy;              % scale frequencies
h     = bar(yi,yfi,'FaceColor',colorVar,'EdgeColor','none');

% Customize plot:
title('Normalized Histogram');
xTxt = sprintf('Values  (n = %d)',size(Y,1));
xlabel(xTxt);
ylabel('Relative Frequency');
% --------------------------
