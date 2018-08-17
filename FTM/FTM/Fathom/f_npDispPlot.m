function hdl = f_npDispPlot(result,conf,iter,hull,peel)
% - plot distance-based measures of Homogeneity in Multivariate Dispersion
%
% USAGE: f_npDispPlot(result,conf,iter,hull,peel);
%
% result  = structure of results obtained from f_npDispPlot
% conf    = confidence interval enclosed by ellipse  (default = 0.95)
% iter    = # iterations for bootstrap resampling    (default = 0)
% hull    = create convex hull                       (default = 0)
% peel    = peeling level for convex hull            (default = 0)
%
% hdl     = handle to axis of group plot for customizing ledgend
%            e.g., legend(hdl,cellArray);
%
% SEE ALSO: f_npDisp

% -----Notes:-----
% If CONF = 0, confidence ellipses are not created.
%
% If ITER = 0, parametric confidence ellipses are generated for each group,
% otherwise a boostrapped version is computed. ITER is recommended to be at
% least 1000 for the latter.

% -----References:-----
% modified after code for f_cdaPlot

% -----Author:-----
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2010: updated documentation, replaced 'unique' with 'f_unique'

% -----Unwrap variables from structure:-----
scores    = result.scores;
centroids = result.C;
grps      = result.grps;
gLabels   = result.gLabels;
% ------------------------------------------

% -----Set defaults and check input:-----
if (nargin < 2), conf   = 0.95; end; % default 95% confidence ellipse
if (nargin < 3), iter   = 0;    end; % default 0 bootstrap samples
if (nargin < 4), hull   = 0;    end; % default no convex hull
if (nargin < 5), peel   = 0;    end; % default no peeling
% -----------------------------------------------

% -----Create Plot:-----
figure('Name','Homogeneity of Multivariate Dispersion');
set(gcf,'color','w');    % set bg color to white
hold on;

uGrps   = f_unique(grps); % unique groups, unsorted
noGrps  = length(uGrps);  % # unique groups   

hdl = zeros(1,noGrps); % preallocate
% plot points
for j = 1:noGrps
   gRows = find(grps==uGrps(j)); % get row indices for each group
   hdl(j) = plot(scores(gRows,1),scores(gRows,2),f_symb(j),'MarkerFaceColor',f_rgb(j),'MarkerEdgeColor',f_rgb(j));
end

% Create legend:
legend(hdl,gLabels);

% plot centroids afterwards, so won't be behind any points:
for j = 1:noGrps
   h = text(centroids(j,1),centroids(j,2),gLabels(j));
   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',10,...
      'Color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',0.5);
end      

xText = ['PCoA Axis I'];
yText = ['PCoA Axis II'];   
xlabel(xText);
ylabel(yText);   

% Plot confidence ellipses:
if (conf>0)
   if (iter>0)
      fprintf('\nBootstrapping the data %d times...\n',iter-1);
   end
   for j = 1:noGrps
      gRows = find(grps==uGrps(j)); % get row indices for each group
      f_confEllipse(scores(gRows,1:2),conf,iter,1,f_rgb(j));
   end
end

% Plot convex hulls:
if (hull>0)
   for j = 1:noGrps
      gRows = find(grps==uGrps(j)); % get row indices for each group
      f_convHull(scores(gRows,1:2),peel,1,f_rgb(j));
   end
end

% Adjust figure properties:
axis tight;
axis(1.05*(axis)); % increase bounds of axis
axis equal
box on;
f_origin('hv'); % mark origin
