function hdl = f_cdaPlot(result,conf,iter,hull,peel,vec,offset,labels,fmt)
% - plot results from a Canonical Discriminant Analysis
%
% USAGE: hdl = f_cdaPlot(result,conf,iter,hull,peel,vec,offset,labels,fmt);
%
% result = structure of results obtained from f_cda
% conf   = confidence interval enclosed by ellipse              (default = 0.95)
% iter   = # iterations for bootstrap resampling                   (default = 0)
% hull   = create convex hull                                      (default = 0)
% peel   = peeling level for convex hull                           (default = 0)
% vec    = include biplot vectors                                  (default = 0)
%          use vec = 2 to create a separate figure
% offset = biplot vector label offset factor ranging from 0-1      (default = 0)
% labels = cell array of variable labels                  (if empty, autocreate)
%           e.g., sLabels = {'sal' 'tmp' 'elev'};
% fmt    = format of labels; 'tex' = TeX formatting (default) or 'none'
% 
% hdl    = handle to axis of group plot for customizing ledgend
%           e.g., legend(hdl,cellArray);
%
% SEE ALSO: f_cda, f_rdaPlot

% -----Notes:-----
% If CONF = 0, confidence ellipses are not created.
%
% If ITER = 0, parametric confidence ellipses are generated for each group,
% otherwise a boostrapped version is computed. ITER is recommended to be at
% least 1000 for the latter.
%
% Vectors indicating variable importance are only created if METHOD 2
% was used in F_CDA (i.e., discriminant vs. classification functions were
% created)
%
% For convex hulls, PEEL must range from 0 to 3 (see f_convHull)

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, Dec-2003 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jan-2004: vectors now optional with scale = 0;
% Feb-2004: added optional convex hulls
% Aug-2004: set bg color to white
% Jul-2004: changed num2Str to num2str
% Sep-2006: edited legend to work with ver. 2006a; misc edits recommended
%           by MLINT
% Oct-2006: fixed problem with legend
% Jun-2007: no longer places title on graph (put this in your figure caption);
%           added support to allow customizing ledgend after function is run
% Oct-2010: replaced unique with f_unique
% Jan-2011: added new method for applying OFFSET; support for FMT
% Mar-2011: trim Cvects to 2 dimensions

% -----Unwrap variables from structure:-----
scores    = result.scores;
Cvects    = result.Cvects(:,1:2);
centroids = result.centroids;
amongVar  = result.amongVar;
grps      = result.y;
gLabels   = result.gLabels;
method    = result.method;
% ------------------------------------------

% -----Set defaults and check input:-----
if (nargin < 2), conf   = 0.95; end; % default 95% confidence ellipse
if (nargin < 3), iter   = 0;    end; % default 0 bootstrap samples
if (nargin < 4), hull   = 0;    end; % default no convex hull
if (nargin < 5), peel   = 0;    end; % default no peeling
if (nargin < 6), vec    = 0;    end; % default no biplot vectors
if (nargin < 7), offset = 0;    end; % set default offset to 0
if (nargin < 8), labels = cellstr(num2str((1:size(Cvects,1))'))'; end; % default Y labels
if (nargin < 9), fmt    = 'tex'; end % default use TeX formatting

% Check offset:
if (offset<0) || (offset>1)
   error('Offset must be in the range 0-1!')
end

% If labels are not cell arrays, try forcing them:
if iscell(labels)<1, labels = num2cell(labels); end;

labels = labels(:); % force as row vector
if (size(Cvects,1) == size(labels,1))<1
   error('The # of rows in RESULT.CVECTS & LABELS must be equal!');
end
% -----------------------------------------------

% -----Scale Cvects, Scores, & Centroids:-----
% - want all coordinates to fit within bounds of Cvects (after eq. 1.10 in L&L, 1998):
% Cvects:
maxVar = max(abs(Cvects(:)));
Cvects = (Cvects./repmat(max([max(abs(Cvects(:,1)));max(abs(Cvects(:,2)))]),...
   size(Cvects)))*maxVar;

% Scores & Centroids:
fac       = max([max(abs(scores(:,1)));max(abs(scores(:,2)))]);
scores    = (scores./repmat(fac,size(scores)))*maxVar;
centroids = (centroids./repmat(fac,size(centroids)))*maxVar;
% --------------------------------------------

% -----Create Canonical Plot:-----
Hm = figure('Name','Canonical Discriminant Analysis');
set(gcf,'color','w');    % set bg color to white
hold on;

uGrps  = f_unique(grps); % unique groups, unsorted
noGrps = length(uGrps);  % # unique groups   

% title('\bfCanonical Discriminant Analysis');
hdl = zeros(1,noGrps); % preallocate
% plot points
for j = 1:noGrps
   gRows  = find(grps==uGrps(j)); % get row indices for each group
   hdl(j) = plot(scores(gRows,1),scores(gRows,2),f_symb(j),'MarkerFaceColor',f_rgb(j),'MarkerEdgeColor',f_rgb(j));
end

% Create legend:
legend(hdl,gLabels);

% Plot centroids afterwards, so won't be behind any points:
for j = 1:noGrps
   h = text(centroids(j,1),centroids(j,2),gLabels(j));
   set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',10,...
      'Color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',0.5);
end      

% Percent of total variation among GROUPS accounted for
% by Axis 1 & 2:
can1 = sprintf('%4.2f',(amongVar(1))*100);
can2 = sprintf('%4.2f',(amongVar(2))*100);

xText = ['Canonical Axis I  (' num2str(can1) ' %)'];
yText = ['Canonical Axis II (' num2str(can2) ' %)'];   
xlabel(xText);
ylabel(yText);   

% ----- Plot confidence ellipses:-----
if (conf>0)
   if (iter>0)
      fprintf('\nBootstrapping the data %d times...\n',iter-1);
   end
   for j = 1:noGrps
      gRows = find(grps==uGrps(j)); % get row indices for each group
      f_confEllipse(scores(gRows,1:2),conf,iter,1,f_rgb(j));
   end
end

% ----- Plot convex hulls: -----
if (hull>0)
   for j = 1:noGrps
      gRows = find(grps==uGrps(j)); % get row indices for each group
      f_convHull(scores(gRows,1:2),peel,1,f_rgb(j));
   end
end


% ----- Biplot vectors (only for method 2 of f_cda):-----
if (method==2) && (vec>0)
   n = size(Cvects,1); % get # of variables
         
   % Get length of each vector:
   eDis = f_dis([0 0;Cvects],'euc');
   eDis = eDis(2:end,1);
   
   % Set delta inversely proportional to length of longest vector:
   ratio  = (abs(eDis./ repmat(max(eDis),size(eDis))));
   ratio  = ones(size(ratio))./ratio;
   offset = repmat(offset,size(eDis));
   delta  = 1 + (offset .* ratio);
   
   % Optionally create a separate figure:
   if (vec>1)
      Hv = figure('Name','Biplot Vectors');
      set(gcf,'color','w');    % set bg color to white
      hold on;
      xText = ['Canonical Axis I  (' num2str(can1) ' %)'];
      yText = ['Canonical Axis II (' num2str(can2) ' %)'];
      xlabel(xText);
      ylabel(yText);
   end
   
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[Cvects(j,1) Cvects(j,2)],'size',[0.125]*0.75,'angle',20,'Color','b')
      
      % Label vectors:
      h = text(Cvects(j,1)*delta(j),Cvects(j,2)*delta(j),labels(j));
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','b','Interpreter',fmt);
   end
   
   % Set figure properties of vector plot:
   if (vec>1)
      figure(Hv); 
      sub_figProp;
   end
end

% Set figure properties of main plot:
figure(Hm); % switch back to main figure
sub_figProp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SUBFUNCTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sub_figProp
% - subfunction to set custom figure properties
axis tight;
axis(1.1*(axis)); % increase bounds of axis
axis equal
box on;
f_origin('hv'); % mark origin
f_arrow;        % re-draw arrows after changing axis scaling











