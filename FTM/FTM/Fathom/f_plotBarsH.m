function f_plotBarsH(x,grp,rel,labels,err,width)
% - horizontal bar graphs
%
% USAGE: f_plotBarsH(x,grp,rel,labels,err,width)
%
% x       = abundance data (columnwise)
% grp     = col vector of integers specifying group membership
% rel     = use relative abundance instead of group mean (default = 0)
% labels  = cell array of variable labels            (if empty, autocreate)
%            e.g., sLabels = {'sal' 'tmp' 'elev'};
% err     = add error bars, where 0 = none (default), 1 = SD, 2 = SE
% width   = width of  bars (0-1)
%
% SEE ALSO: f_plotBarsV
%

% -----Author:-----
% by David L. Jones, 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% - after f_plotHBars/f_plotBarsV by D. Jones, June-2007

% -----Check input & set defaults:-----
if (nargin <3),rel    = 0; end; % default no relative abundance
if (nargin <4)                  % default Y labels
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


% Create Bar graph:
barh(flipud((1:n)'),res,width);
set(gca,'FontSize', 8,'TickDir','in','YTickLabel',labels);

% Set Bar Colors:
if (nt<2)
   colormap(repmat(((0:nt)')/nt,1,3));
else
   % For grouped bar graphs, 1st is white, 2nd is black:
   colormap(repmat(((nt:-1:0)')/nt,1,3));
end

% Add error bars:
if (err>0) && (rel<1) % no error for relative abundance
   if (err == 1) % standard deviation
      f_plotError(flipud((1:n)'),res,gStdv,2);
   elseif (err == 2) % standard error
      if nt==1     % 1 treatment
         f_plotError(flipud((1:n)'),res,gSE,2);
      elseif nt==2 % 2 treatments
         f_plotError(flipud(((1:n)')-0.125),res(:,1),gSE(:,1),2);
         f_plotError(flipud(((1:n)')+0.125),res(:,2),gSE(:,2),2);
      end
   end
end
% reverse stacking order so error bars are behind:
set(gca,'children',flipud(get(gca,'children')));
%
% bring tick marks & grid line to the front:
set(gca,'Layer','top');






