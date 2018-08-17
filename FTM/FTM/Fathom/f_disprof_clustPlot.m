function f_disprof_clustPlot(result,txt,top,scale)
% - plot dendrogram of a DISPROF-based custer analysis
%
% USAGE: f_disprof_clustPlot(result,txt,top,scale);
%
% result = structure of results obtained from 'f_disprof_clust' function
% txt    = cell array of row labels, if empty autocreate
%          e.g., txt = {'A' 'B' 'C'};
% top    = orientation is top-to-bottom vs. right-to-left          (default = 1)
% scale  = apply scaling factor to size of axis and tick labels    (default = 1) 
%
% SEE ALSO: f_disprof_clust

% -----Author:-----
% by David L. Jones, Jan-2014
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% Apr-2014: added option for TOP and SCALE; tick marks are plotted on one axis

% -----Set defaults & check input:-----
if (nargin < 2), txt   = []; end % default create grouping vector
if (nargin < 3), top   = 1;  end % default top-to-bottom orientation
if (nargin < 4), scale = 1;  end % default apply scaling factor of 1 to fonts

% Set default grouping vector:
if isempty(txt)
   txt = f_num2cell([1:numel(result.grp)]);
end

% Set orientation of dendrogram:
if (top<1)
   orientVar = 'right';
else
   orientVar = 'top';
end

% Check scaling factor:
if (scale<=0), error('SCALE must be greater than 0!'); end
% -------------------------------------

% Make sure labels are of compatible size:
txt = txt(:); % force cell array into a column

if numel(txt) ~= size(result.grp,1)
   error('TXT & RESULT.Z don''t have the same # of rows!')
end

% Create dendrogram:
figure;
set(gcf,'color','w'); % set background color to white
H = dendrogram(result.Z,0,'LABELS',txt,'ORIENTATION',orientVar);

% Flip the order of graphics handles:
H = flipud(H);

% Set the color:
set(H(result.colS==1),'Color','k');          % heterogeneous clusters
set(H(result.colS==0),'Color',[1 1 1]*0.75); % homogeneous clusters

% Increase/decrease font size (applies to axis labels and tic labels):
fontVar = get(gca,'FontSize');
set(gca,'FontSize',fontVar * scale)

% Set tic direction:
set(gca,'TickDir', 'out');

% Reduce tick length by half:
ticLen = get(gca,'TickLength');
set(gca,'TickLength',[ticLen(1)/2 ticLen(2)]);

% Customize plot:
title('\bfDISPROF-based Cluster Analysis');
if (top<1)
   % Set X-axis label:
   xlabel('\bfDissimilarity');
   
   % Add some padding:
   xtickVar = get(gca,'xtick');
   xlimVar  = xlim;
   xlim([xlimVar(1) - 0.1*mean(diff(xtickVar)) xlimVar(2)]);

else
   % Set Y-axis label:
   ylabel('\bfDissimilarity');
   
   % Add some padding:
   ytickVar = get(gca,'ytick');
   ylimVar  = ylim;
   ylim([ylimVar(1) - 0.1*mean(diff(ytickVar)) ylimVar(2)]);
end

% -----Customize tic marks so they're just on 1 side:-----
% http://stackoverflow.com/questions/15553720/matlab-remove-only-top-and-right-ticks-with-leaving-box-on
a = gca; % get handle to current axes
set(a,'box','off','color','none'); % set box property off, remove bg color
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set original axes as active
axes(a)
% link axes in case of zooming
linkaxes([a b])
