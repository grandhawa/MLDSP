function f_grpPlot(scores,grp,gLabels,sm)
% - group plotting function
% 
% USAGE: f_grpPlot(scores,grp,gLabels,sm)
% 
% scores  = coordinates
% grp     = column vector of whole numbers indicting group membership
% gLabels = cell array of group labels; if empty, autocreate
%            e.g., gLabels = {'one' 'two' 'three'};
% sm       = use spatial median instead of centroid                  (default = 0)

% -----Author:-----
% by David L. Jones, Dec-2012
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), gLabels  = []; end % default autogenerate group labels
if (nargin < 4), sm   = 0;      end % default use centroid vs. spatial median

% Check size of input:
n = size(scores); % # obs
if n ~= size(grp,1), error('SCORES & GRP need same # of rows'); end

% Autogenerate group labels:
if isempty(gLabels)
   gLabels = f_num2cell(f_unique(grp));
end

% If labels are not cell arrays, try forcing them:
if iscell(gLabels)<1, gLabels = num2cell(gLabels); end

% Set parameters:
fmt = 'none';
% -------------------------------------

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
scores = (scores./repmat(max([max(abs(scores(:,1)));max(abs(scores(:,2)))]),...
   size(scores)));

% Find group centroids:
if (sm>0)
   [centroids,SE] = f_centroid(scores,grp,2); % spatial median
else
   [centroids,SE] = f_centroid(scores,grp,1); % centroid
end

% -----Create Plot:-----
figure;
set(gcf,'color','w');    % set bg color to white
hold on;

uGrp = f_unique(grp); % unique groups, unsorted
nGrp = length(uGrp);  % # unique groups

% Plot sites:
hdl = zeros(1,nGrp); % preallocate
for j = 1:nGrp
   gRows  = grp==uGrp(j); % get row indices for each group
   hdl(j) = plot(scores(gRows,1),scores(gRows,2),'o','MarkerFaceColor',...
      [1 1 1]*1,'MarkerEdgeColor','none');
end

% Plot standard error bars:
f_plotError(centroids(:,1),centroids(:,2),SE(:,1),1); % Axis I
f_plotError(centroids(:,1),centroids(:,2),SE(:,2),2); % Axis II

% Plot centroids:
for j = 1:nGrp
   h = text(centroids(j,1),centroids(j,2),gLabels(j));
   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',10,...
      'Color',[0 0 0],'BackgroundColor','none','Interpreter',fmt);
%    set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',10,...
%       'Color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',0.25,'Interpreter',fmt);
end

% Customize plot:
axis tight;
axis(1.05*(axis)); % increase bounds of axis
axis equal;
box on;
f_origin('hv');


