function [h,idx] = f_convHull(xy,peel,plt,cSpec)
% - convex hull for bivariate data
%
% USAGE: f_convHull(xy,peel,plt,cSpec)
%
% xy    = 2 column matrix of input data
% peel  = peeling level                           (default = 0)
% plt   = add ellipse to current plot             (default = 0)
% cSpec = line spec defining ellipse color        (default = 'r');
% 
% h   = handle to plotted convex hull
% idx = index of vectors of convex hull
%
% % SEE ALSO: f_confEllipse

% -----Notes:-----
% This function is used to create and plot convex hulls, with
% optional peeling. Convex hulls are useful in exploratory data analysis to
% help distinguish groups, outliers, and general shapes in the data. Peeled
% convex hulls are equivalent to bivariate smoothers. Unlike confidence
% ellipses, no assumptions concerning the underlying distribution of the
% data are made.

% -----Author:-----
% by David Jones, Feb-2004
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 2), plt   = 0;   end; % default no peeling
if (nargin < 3), plt   = 0;   end; % default don't plt ellipse
if (nargin < 4), cSpec = 'r'; end; % default color 'r'


if (size(xy,2) ~= 2)
   error('XY must be a 2 column matrix')
end

if (peel < 0) || (peel > 3)
   error('PEEL must be between 0 and 3');
end
% ---------------------------------------

% Parse input data:
x = xy(:,1);
y = xy(:,2);

% Get index of vectors of convex hull:
idx = convhull(x,y);

% Optionally peel convex hull:
if (peel>0)
   for i=1:peel
      x(idx) = [];            % peel previous hull
      y(idx) = [];            % peel previous hull
      idx    = convhull(x,y); % find new hull
   end
end

% Add convex hull to current plot:
if (plt>0)
   hold on;
   h = line(x(idx),y(idx),'Clipping','off');
   set(h,'Color',cSpec);
else
   h = NaN;
end

