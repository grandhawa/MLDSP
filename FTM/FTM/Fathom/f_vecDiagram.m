function f_vecDiagram(u,v,theta,units,jdate)
% - plot progressive vector diagrams
%
% USAGE: f_vecDiagram(u,v,theta,units,jdate)
%
% u,v    = vector components
% theta  = angle to rotate vectors ccw (default = 0)
% units  = m/s (1) or cm/s (2)         (default = 0)
% jdate  = Julian dates corresponding to u,v
%
% See also: f_vecPlot

% ----- Notes: -----
% This function is used to create progressive vector diagrams from time
% series data of wind or moored current meter velocity vectors. This type
% of diagram is used to produce a Lagrangian display of Eulerian
% measurements.
%
% THETA is an optional angle of rotation that can be applied to rotate
% vectors counter-clockwise (or the coordinate system clockwise) so they, for
% example, align with local isobath coordinate system. Doing so partitions
% the vectors into alongshore (u) and cross-shore (v) components.
%
% UNITS is an optional parameter that allows calculation of the spatial
% units in the plot. A velocity vector specifying 1 m/s covers 3.6 km/hr
% (there is 3600 s in an hour).
%
% JDATE is an optional column vector of Julian dates corresponding to U,V.
% If provided plot symbols will only be drawn for vectors corresponding to
% whole Julian dates.
%
% The aspect ratio of the plot should provide an indication of the relative
% contribution of the U vs. V components to the net displacement. Note that
% the starting position of the diagram is the origin (0,0).

% ----- Author: -----
% by David L. Jones, Dec-2002
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: optional angle of rotation (theta), updated documentation,
%           optionally plot symbols for only whole Julian dates
% Mar-2008: used logical indexing vs. find; set bg color


% ----- Check input & set defaults: -----
u = u(:); % force column vector
v = v(:);

if (size(u,2)~=1)
	error('U & V must be column vectors!');
end

if (size(u) ~= size(v))
	error ('U & V must be same size!');	
end

if (nargin < 3), theta  = 0; end; % no vector rotation
if (nargin < 4), units  = 0; end; % no units provided

if (nargin < 5)
   pltAll = 1;       % plot symbols for all data
else
   pltAll = 0;       % only plot symbols for whole Julian dates
   jdate = jdate(:); % force column vector
   if size(jdate,1) ~= size(u,1)
      error('JDATE must have same # ROWS as U,V!');   
   end
end;

if theta<0, error('THETA must an angle in POSITIVE degrees!'); end
% ---------------------------------------


% -----Optionally rotate vectors counter-clockwise:-----
if theta>0
   % Convert to mag/dir:
   [mag,dir] = f_vecMagDir(u,v);
   
   % Rotate vectors counter-clockwise theta degrees:
   dir = dir + theta;
   dir(dir>360) = dir(dir>360) - 360;
   
   % Convert back to u,v:
   [u,v] = f_vecUV(mag,dir);
end
% ------------------------------------------------------

% Determine spatial scale of diagram:
switch units
case 0 % no units
case 1 % m/s
	u = u*3.6;
	v = v*3.6;
case 2 % cm/s
	u = u*0.036;
	v = v*0.036;
otherwise
	error('Unsupported UNITS');	
end

% nr = size(u,1);

tail = cumsum([[0 0];[u v]]);  % start from the origin
head = [tail(2:end,:);[NaN NaN]];

figure;
set(gcf,'color','w'); % set bg color to white
hold on;

% Plot lines connecting tails to heads:
plot([tail(:,1) head(:,1)]',[tail(:,2) head(:,2)]','b-');

% Plot symbols at tail:
if (pltAll>0) % for all vectors
   plot(tail(:,1),tail(:,2),'r.');
   
else % Only for whole Julian dates
   idx = find(jdate - fix(jdate)==0);
   plot(tail(idx,1),tail(idx,2),'r.');
   
   % Plot Julian dates:
   text(tail(idx,1),tail(idx,2),num2cell(jdate(idx)));
   
end

title('Progressive Vector Diagram');

if (units>0) % spatial
	xlabel('W to E Transport (km)');
	ylabel('S to N Transport (km)');
else % velocity
	xlabel('Velocity of U component');
	ylabel('Velocity of V component');
end

% Adjust appearance of plot:
set(gcf,'color','w');
grid on;
box on;
daspect([1 1 1]);

hold off;