function h = f_plotNeigh(N,plt,colr,txt)
% - plotting function for f_delaunay, f_dnn, f_gabriel, f_mst, & f_relNeigh
%
% USAGE: h = f_plotNeigh(N,plt,colr,txt)
%
% N    = neighbor graph created by f_delaunay, f_dnn, f_gabriel, f_mst, or
%        f_relNeigh
% plt  = plot symbols (= 1, default), labels (= 2), or none (= 0)
% colr = colorSpec for lines indicating neighborhood connections (default = 'k')
% txt  = cell array of plot labels used when plt=2; if empty, autocreate
%          e.g., txt = {'s1' 's2' 's3'};
%
% h = handle to figure axis
%
% SEE ALSO: f_delaunay, f_dnn, f_gabriel, f_mst, f_relNeigh

% -----Author:-----
% by David L. Jones, March-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jan-2013: updated documentation; not supports user-supplied plot labels
%           (txt).

% -----Check inputs and set defaults:-----
if (nargin < 2), plt  = 1;   end % default plot symbols
if (nargin < 3), colr = 'k'; end % colorSpec for
if (nargin < 4), txt  = [];  end % default no plot labels

if isnan(sum(N.dat(:)))
   error('N.dat is empty, re-create with original coordinates!')
end

if ((size(N.dat,1) ~= size(N.dis,2)))
   error('Size mismatch between N.DIS & N.DAT!')
end

if plt==2
   % Autocreate yLabels:
   if isempty(txt)
      txt = cellstr(num2str([1:size(N.dat,1)]')); % create some labels:
   end
         
   % If labels are not cell arrays, try forcing them:
   if iscell(txt)<1, txt = num2cell(txt); end
   
   % Make sure labels are of compatible size:
   txt = txt(:);
   if size(txt,1) ~= size(N.dat,1)
      error('Size mismatch between N.dat and txt!')
   end
end
% ----------------------------------------

% Extract components:
crds = N.dat;
con  = N.con;

% Set plot options:
set(gcf,'color','w'); % set bg color to white
hold on;
box on;

% Plot lines indicating connecting neighbors:
plot([crds(con(:,1),1) crds(con(:,2),1)].',[crds(con(:,1),2) crds(con(:,2),2)].',...
   '-', 'Color', colr);

if (plt==1)     % Overlay symbols:
h = plot(crds(:,1),crds(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor',...
      'none','MarkerSize',4);

elseif (plt==2) % Overlay text labels:
   h   = text(crds(:,1),crds(:,2),txt);
   set(h,'HorizontalAlignment','center','Color','k','FontSize',12,...
      'FontWeight','bold','Interpreter','none')

elseif (plt==0)
   % Do nothing

else
   error('PLT must be 0, 1, or 2')
end

% More plot options:
daspect([1 1 1]);
