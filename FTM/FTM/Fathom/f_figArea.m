function [area,h2] = f_figArea(h);
% - finds the area of figure objects
%
% USAGE: [area,h2] = f_figArea(h);
%
% h    = handle of figure object(s)
% area = area for object specified by h
% h2   = figure handles sorted by area (ascending)

% -----Notes:-----
% This function is used to compute the Area of figure objects
% in a plot and optionally returns a list of handles sorted
% ascending by Area. This is particularly useful when you've
% plotted a number of patches but, because of the stacking order,
% the smaller ones are obscured by the larger.

% -----Examples:-----
% h = get(gca,'Children');  % get stacking order of handles
% [null,h2] = f_figArea(h); % sort by area
% set(gca,'Children',h2);   % re-stack

% -----References:-----
% formula for area after comments in mu_coast.m
% from M_Map toolbox by Rich Pawlowicz's <rich@ocgy.ubc.ca>
% available from:
% http://www2.ocgy.ubc.ca/~rich/map.html

% -----Author(s): -----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

n = size(h,1); % get # of objects

% Find the Area of each object:
for i = 1:n
   x = get(h(i),'xdata'); x = x(:);
   y = get(h(i),'ydata'); y = y(:);
   nl = length(x);
   area(i) = sum(diff(x).*(y(1:nl-1) + y(2:nl))/2);
end

area = area(:); % column vector

res = sortrows([h area],2);  % sort handles by area
h2  = (res(:,1));            % sorted handles





