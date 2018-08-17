function f_origin(where,type,level)
% - mark the origin of a graph with horizontal/vertical line(s)
%
% USAGE: f_origin('where','type',level)
%
% where = 'h' (horizontal), 'v' (vertical), 'hv' (both)
% type  = linespec (default = ':')
% level = gray level, smaller values are darker (default = 0.70)

% -----Notes:-----
% This function is used to mark the origin of the current graph with a
% horizontal and/or vertical (dotted) line.

% -----Author:-----
% by Dave Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jan-2009: added support for user-specified gray level

% -----Check input & set defaults:-----
if (nargin < 2), type  = ':'; end  % dotted line by default
if (nargin < 3), level = 0.70; end % default gray level

where = lower(where); % force small caps

% allow either designation
if (where=='vh')
   where = 'hv';
end

al = axis;                % get current axis limits
h  = get(gca,'Children'); % get original stacking order of handles

% Create lines:
switch where
   case 'h' % horizontal line:
      h_new = line([al(1);al(2)],[0;0],'LineStyle',type,'Color',[1 1 1]*level);
      
   case 'v' % vertical line:
      h_new = line([0;0],[al(3);al(4)],'LineStyle',type,'Color',[1 1 1]*level);      
           
   case 'hv' % both:
      h_h = line([al(1);al(2)],[0;0],'LineStyle',type,'Color',[1 1 1]*level);
      h_v = line([0;0],[al(3);al(4)],'LineStyle',type,'Color',[1 1 1]*level); 
      h_new = [h_v;h_h];
             
   otherwise
      error('WHERE must be ''h'', ''v'', or ''hv'' !')
end

% Adjust stacking order (new lines are on botom):
h2 = [h;h_new];
set(gca,'Children',h2);

% Force previous axis limits:
axis(al);

