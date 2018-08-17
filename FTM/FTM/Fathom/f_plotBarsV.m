function h = f_plotBarsV(x,grp,rel,labels,err,width)
% - vertical bar graphs
%
% USAGE: h = f_plotBarsV(x,grp,rel,labels,err,width);
%
% x       = abundance data (columnwise)
% grp     = col vector of integers specifying group membership
% rel     = use relative abundance instead of group mean (default = 0)
% labels  = cell array of variable labels            (if empty, autocreate)
%            e.g., sLabels = {'sal' 'tmp' 'elev'};
% err     = add error bars, where 0 = none (default), 1 = SD, 2 = SE
% width   = width of  bars (0-1)
% 
% h       = handle to bar graphics
%
% SEE ALSO: f_plotBarsH
%

% -----Author:-----
% by David L. Jones
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% - after f_plotVBars by D. Jones, April-2003

% Jun-2007: modified colormap to accept > 2 groups, added support for SE,
%           removed FIGURE and SET commmands so will work inside a SUBPLOT,
%           added support for error bars when the number of treatments = 2.
% 
% Oct-2009: updated to simplier way to assign colormap, return h

% -----Notes:-----
% Recommend width = 0.25 for 1 factor plots and 0.7 for 2 factor plots.

% -----Check input & set defaults:-----
if (nargin <3),rel    = 0; end; % default no relative abundance
if (nargin <4)                  % default X labels
   labels = cellstr(num2str((1:size(unique(grp),1))'))'; end
if (nargin <5),err = 0; end      % default no error bars
if (nargin <6),width = 0.5; end  % default width of bars
if (width>1), error('WIDTH is > 1!'); end
% -------------------------------------

if (rel<1)
   [res,gStdv,gSE] = f_grpMean(x,grp);
else
   res = f_grpRel(x,grp);
end

n  = size(res,1);   % # of groups
nt = size(x,2);     % # of treatments

% figure;               % open new figure window
% set(gcf,'color','w'); % set bg color so printed/exported lakes will be white

% Create Bar graph:
h = bar((1:n)',res,width); % <== use for 1 Factor plots

set(gca,'FontSize', 8,'TickDir','out','XTickLabel',labels);


% % Set Bar Colors:
% if (nt<2)
%    colormap(repmat(((0:nt)')/nt,1,3));
% else
%    % For grouped bar graphs, 1st is white, 2nd is black:
%    colormap(repmat(((nt:-1:0)')/nt,1,3));
% end

% Use gray-scale colormap:
colormap(f_colormap(nt));

% Add error bars:
if (err>0) && (rel<1) % no error for relative abundance
   if (err == 1) % standard deviation
      f_plotError((1:n)',res,gStdv);
   elseif (err == 2) % standard error
      if nt==1     % 1 treatment
         f_plotError((1:n)',res,gSE,1,1.5);                 % width = 1.5
      elseif nt==2 % 2 treatments
         f_plotError(((1:n)')-0.125,res(:,1),gSE(:,1),1,2); % width = 1.5
         f_plotError(((1:n)')+0.125,res(:,2),gSE(:,2),1,2); % width = 1.5
      end
   end
end
% Reverse stacking order so error bars are behind:
set(gca,'children',flipud(get(gca,'children')));
%
% Bring tick marks & grid line to the front:
set(gca,'Layer','top');






