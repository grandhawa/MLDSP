function [hdl,yCor,centroids,crds] = f_capPlot(result,gLabels,sLabels,Y,yLabels,offset,fmt,wascores,ctr)
% - plot results from a CAP analysis
%
% USAGE: [hdl,yCor,centroids,crds] = f_capPlot(result,gLabels,sLabels,Y,yLabels,offset,fmt,wascores,ctr);
%
% result   = structure of results obtained from f_cap
% gLabels  = cell array of group labels; if empty, autocreate
%            e.g., gLabels = {'one' 'two' 'three'};
%
% sLabels  = cell array of site labels; if empty, plot as colored circles
%            e.g., sLabels = {'A' 'A' 'A' 'B' 'B' 'C'};
%
% Y        = matrix of (transformed) response variables
% yLabels  = cell array of Y labels; if empty, autocreate
%           e.g., yLabels = {'sp1' 'sp2' 'sp3'};
%
% offset   = yLabel offset factor ranging from 0-1                 (default = 0)
% fmt      = format of labels; 'tex' = TeX formatting (default) or 'none'
% wascores = plot Y's weighted average scores                      (default = 0)
% ctr      = plot centroids                                        (default = 0)
%
% hdl      = handle to axis of group plot for customizing ledgend
%            e.g., legend(hdl,cellArray);
% yCor      = biplot vector coordinates
% centroids = centroids scaled to fit in -1 to 1 bounds
% crds      = siteScores scaled to fit in -1 to 1 bounds
%
% SEE ALSO: f_cap, f_cdaPlot, f_rdaPlot

% -----Author:-----
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% Jul-2011: updated documentation
% Nov-2011: updated documentation
% Apr-2012: don't create a legend if sLabels are provided; updated
%           documentation
% May-2012: don't title the CAP plot
% Oct-2012: now optionally plots wascores
% Nov-2012: removed error when plotting jittered centroids; plotting
%           centroids is now optional
% Mar-2013: corrected misspelling in error message
% Nov-2013: updated documentation
% Mar-2014: now outputs scaled siteScores for use outside the function

% -----Unwrap variables from structure:-----
siteScores = result.siteScores;
centroids  = result.centroids;
canVar     = result.canVar;
grp        = result.grp;
% ------------------------------------------

% -----Set defaults and check input:-----
if (nargin < 2), gLabels  = [];    end % default autogenerate group labels
if (nargin < 3), sLabels  = [];    end % default don't plot site labels
if (nargin < 4), Y        = [];    end % no biplot
if (nargin < 5), yLabels  = [];    end % default Y labels
if (nargin < 6), offset   = 0;     end % set default offset to 0
if (nargin < 7), fmt      = 'tex'; end % default use TeX formatting
if (nargin < 8), wascores = 0;     end % default don't plot wascores
if (nargin < 9), ctr      = 0;     end % default don't plot centroids

% Check offset:
if (offset<0) || (offset>1)
   error('Offset must be in the range 0-1!')
end

% Autogenerate group labels:
if isempty(gLabels)
   gLabels = cellstr(num2str(f_unique(grp))); % create group labels for f_cdaPlot
end

% If labels are not cell arrays, try forcing them:
if iscell(gLabels)<1, gLabels = num2cell(gLabels); end

% Create/check yLabels:
if ~isempty(Y)
   % Get # response variables:
   ncY = size(Y,2);
   
   % Autocreate yLabels:
   if isempty(yLabels)
      yLabels = cellstr(num2str((1:ncY)'))';
   end
   
   % If labels are not cell arrays, try forcing them:
   if iscell(yLabels)<1, yLabels = num2cell(yLabels); end
   
   % Make sure labels are of compatible size:
   yLabels = yLabels(:);
   if size(yLabels,1) ~= ncY
      error('Size mismatch of yLabels and Y!')
   end
end
% ---------------------------------------

% If only 1 canonical axis, add 2nd jittered axis:
if (size(siteScores,2)==1)
   a          = min(siteScores) * 0.25;
   b          = max(siteScores) * 0.25;
   c          = min(centroids)  * 0.25;
   d          = max(centroids)  * 0.25;
   jit        = a + (b-a).*rand(size(siteScores));
   jit_2      = c + (d-c).*rand(size(centroids));
   siteScores = [siteScores jit];
   centroids  = [centroids jit_2];
   canVar     = [canVar;0];
else
   jit = [];
end


% -----Create Canonical Plot:-----
figure('Name','CAP: Canonical Discriminant Analysis');
set(gcf,'color','w');    % set bg color to white
hold on;

uGrp = f_unique(grp); % unique groups, unsorted
nGrp = length(uGrp);  % # unique groups

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
nC         = size(centroids,1);      % get # centroids
siteScores = [centroids;siteScores]; % append centroids
siteScores = (siteScores./repmat(max([max(abs(siteScores(:,1)));max(abs(siteScores(:,2)))]),...
   size(siteScores)));
centroids  = siteScores(1:nC,:); % extract scaled centroids
siteScores(1:nC,:) = [];         % remove them

% Plot sites:
hdl = zeros(1,nGrp); % preallocate
if isempty(sLabels)
   for j = 1:nGrp
      gRows  = grp==uGrp(j); % get row indices for each group
      hdl(j) = plot(siteScores(gRows,1),siteScores(gRows,2),f_symb(j),'MarkerFaceColor',f_rgb(j),'MarkerEdgeColor',f_rgb(j));
   end
   
   % Create legend:
   legend(hdl,gLabels);
   
else
   plot(siteScores(:,1),siteScores(:,2),'.w','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
   text(siteScores(:,1),siteScores(:,2),sLabels,'HorizontalAlignment','center',...
      'VerticalAlignment','middle','Color','b','Interpreter',fmt);
end

% -----Optionally plot centroids:-----
% - plot after siteScores so won't be behind any points:
if (ctr>0)
   for j = 1:nGrp
      h = text(centroids(j,1),centroids(j,2),gLabels(j));
      set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',10,...
         'Color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',0.5,'Interpreter',fmt);
   end
end
% ------------------------------------

% Percent of total variation among GROUPS accounted for
% by Axis 1 & 2:
can1 = sprintf('%4.2f',(canVar(1)/sum(canVar))*100);
can2 = sprintf('%4.2f',(canVar(2)/sum(canVar))*100);

title('\bfCanonical Analysis of Principal Coordinates (CAP)');
xText = ['Canonical Axis I  (' num2str(can1) ' %)'];
if isempty(jit)
   yText = ['Canonical Axis II (' num2str(can2) ' %)'];
else
   yText = 'Jittered Axis';
end

xlabel(xText);
ylabel(yText);

axis tight;
axis(1.05*(axis)); % increase bounds of axis
axis equal;
box on;
f_origin('hv');


% -----Optional species biplot vectors:-----
if (~isempty(Y))
   s = size(siteScores,2);
   % Correlation of original Y variables with each canonical axis:
   yCor = zeros(ncY,s); % preallocate
   for i = 1:s
      for j = 1:ncY
         yCor(j,i) = f_corr(Y(:,j),siteScores(:,i));
      end
   end
   
   % Get length of each vector:
   eDis = f_dis([0 0;yCor(:,1:2)],'euc');
   eDis = eDis(2:end,1);
   
   % Set delta inversely proportional to length of longest vector:
   ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
   ratio = ones(size(ratio))./ratio;
   delta = 1 + (1* repmat(offset,size(eDis)) .* ratio);
   
   
   figure('Name','Correlation Vectors');
   set(gcf,'color','w'); % set bg color to white
   hold on;
   
   nCols  = 2;
   radius = sqrt(2/nCols);
   rectangle('Position',[-radius,-radius,radius*2,radius*2],...
      'Curvature',[1,1],'EdgeColor',[1 1 1]*0.75);
   
   axis([-1 1 -1 1]);
   daspect([1 1 1]);
   
   xlabel('Correlation with Canonical Axis I');
   if isempty(jit)
      ylabel('Correlation with Canonical Axis II');
   else
      ylabel('Jittered Axis');
   end
   
   n = size(yCor,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[yCor(j,1) yCor(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(yCor(j,1)*delta(j),yCor(j,2)*delta(j),yLabels(j));
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter',fmt);
   end
   
   % Customize plot:
   box on;
   axis([-1 1 -1 1]*1.05)
   set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
   f_origin('hv');
   
else
   yCor = NaN;
end
% ------------------------------------------


% -----Optionally plot Y's WA scores (McCune, 1997):-----
if (wascores>0) && ~isempty(Y) % plot
   
   figure('Name','Weighted Averages Scores');
   set(gcf,'color','w'); % set bg color to white
   hold on;   
   
   % Plot sites:
   plot(siteScores(:,1),siteScores(:,2),'o','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
      
   % Percent of total variation among GROUPS accounted for
   % by Axis 1 & 2:
   can1 = sprintf('%4.2f',(canVar(1)/sum(canVar))*100);
   can2 = sprintf('%4.2f',(canVar(2)/sum(canVar))*100);
   
   % Plot WA scores:
   waY = f_wascores(siteScores(:,1:2),Y,0);
   h   = text(waY(:,1),waY(:,2),yLabels);
   set(h,'HorizontalAlignment','center','Color',f_rgb(2),'FontSize',8,...
      'FontWeight','bold', 'Interpreter', fmt);
   
   % Customize plot:
   title('\bfCAP: Weighted Averages Scores');
   xText = ['Canonical Axis I  (' num2str(can1) ' %)'];
   if isempty(jit)
      yText = ['Canonical Axis II (' num2str(can2) ' %)'];
   else
      yText = 'Jittered Axis';
   end
   
   xlabel(xText);
   ylabel(yText);
   
   axis tight;
   axis(1.05*(axis)); % increase bounds of axis
   axis equal;
   box on;
   f_origin('hv');
   
end
% ---------------------------------------------------

% Rename variables for output:
crds = siteScores; % scaled siteScores


