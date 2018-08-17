function f_shadeBox(region,scale)
% - shade subsets of a time series plot
%
% USAGE: f_shadeBox(region,scale) or f_shadeBox('region',scale)
%
% region = 2-d matrix defining regions along the Y-axis
%          to shade (column 1 = start, column 2 = stop)
%
% scale  = scaling factor used in time series plot
%          (default = 1)
%
%  REGION  specifies a variable currently loaded in the Matlab workspace,
%
% 'REGION' specifies a space-delimited ASCII file in the Matlab path.
%
% See also: f_vecPlot

% ----- Notes: -----
% This function is used to shade subsets of a time series plot in order to
% highlight specific time periods. It should be called after creating a time
% series plot and should use the same scaling factor.

% ----- Author: -----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: updated documentation
% Jan-2010: changed '|' to '||'

% ----- Set defaults & check input: -----
if (nargin < 2), scale = 1; end; % no scaling by default

if (scale==0)
	error('You cannot scale by 0');
end
% ---------------------------------------

% ----- Load data & check input: -----
if (ischar(region)) % region specifies data file
	region = textread(region,'','delimiter',' ');
end

[nr,nc] = size(region);

if (nc<2) || (nc>2)
	error('REGION must be a 2-d matrix')
end
% -----------------------------------

% Scale:
region = region/scale;

% Define color variable:
colorVar = [1 1 1]*0.90; % 10% gray

% Get Y-axis bounds:
theAxis = axis;
yMin    = theAxis(3);
yMax    = theAxis(4);

% Define coordinates of box to shade:
for i=1:nr
	shade{i} = [region(i,1) yMin; region(i,1) yMax; region(i,2) yMax; region(i,2) yMin];
end

% Plot filled rectangles as patches:
hold on;
for i=1:nr
	h = fill(shade{i}(:,1),shade{i}(:,2),colorVar);
	set(h,'EdgeColor','none');
end

% Reverse stacking order so patches are behind:
set(gca,'children',flipud(get(gca,'children')));

% Bring tick marks & grid line to the front:
set(gca,'Layer','top');

hold off;

