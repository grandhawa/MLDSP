function [hdl,biplotX,scores] = f_RFvis(model,conf,iter,offset,fmt,idx)
% - visualize a Random Forest
%
% USAGE: [hdl,scores] = f_RFvis(model,conf,iter,offset,fmt,idx);
%
% model  = structure of results obtained from f_RFclass
% conf   = confidence interval enclosed by ellipse               (default = 0.95)
% iter   = # iterations for bootstrap resampling                    (default = 0)
% offset = label offset factor ranging from 0-1                    (default = 0)
% fmt    = format of labels; 'tex' = TeX formatting (default) or 'none'
% idx    = index to biplot vectors to plot                    (default plot all)
%
% hdl     = handle to axis of group plot for customizing ledgend
%           e.g., legend(hdl,cellArray);
% biplotX = coordinates of biplot vectors
% scores  = coordinates of plotted scores for use in other routines
%
% SEE ALSO: f_RFclass, f_RFclassPredict, f_RFplot

% -----Notes:-----
% For UNSUPERVISED LEARING: there is no 'a priori' class structure imposed on
% the data and axes are labelled according to the proportion of the total
% variation in the data they account for. This is analagous to an
% unconstrained ordination such as PCoA, though seems to pick out latent
% structue in the data much better than a tradition PCoA or PCA.
%
% For CANONIAL ANALYSIS: an 'a priori' class structure is imposed on the
% data and axes are labelled according to the proprotion of the
% within-group variation in the data explained by each axis. This is
% analagous to a constrained ordination, such as db-RDA.

% -----References:-----
% Shi, T. and S. Horvath. 2006. Unsupervised learning with Random Forest
% predictors. Journal of Computational and Graphical Statistics 15(1):
% 118-138.

% -----Author:-----
% by David L. Jones, Aug-2009
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Oct-2010: unique replaced by f_unique; support for biplot vectors
% Nov-2010: constrained PCoA now performed via db-RDA which properly trims
%           canonical axes, correctly calculates proportion of among-group
%           variation explained by each axis, etc.
% Aug-2011: PCoA scores and biplot vectors are automatically scaled; biplot
%           vector created as separate figure; now returns biplot
%           coordinates
% Feb-2015: edited to work with changes in f_pcoa

% -----Set defaults and check input:-----
if (nargin < 2), conf   = 0.95;  end % default 95% confidence ellipse
if (nargin < 3), iter   = 0;     end % default 0 bootstrap samples
if (nargin < 4), offset = 0;     end % default is no label offset
if (nargin < 5), fmt    = 'tex'; end % default use TeX formatting
if (nargin < 6), idx    = [];    end % default initialize as empty

% Check offset:
if (offset<0) || (offset>1)
   error('Offset must be in the range 0-1!')
end

% Check for proximities:
if isnan(model.prox)
   error('No proximities exist, re-run f_RFclass with sim=1!')
end
% ---------------------------------------

% Extract data from RANDOM FOREST:
dis    = sqrt(1 - model.prox); % RF dissimilarities (Shi & Horvath, 2006)
nClass = model.nClass;
grps   = model.Y;
X      = model.X;     % transformed predictors
X_txt  = model.X_txt; % cell array of variable names
Y_txt  = model.Y_txt;
clear model;

% Plot all biplot vectors by default
if (isempty(idx))
   idx = 1:size(X,2);
end


% Setup groups:
uGrps   = unique(grps);  % unique groups, sorted
noGrps  = length(uGrps); % # unique groups

if noGrps==1 % Unsupervised learning (unconstrained PCoA):
   pcoa     = f_pcoa(dis,0,1,0); % scale eigenvectors, no negative eigenvalues
   scores   = pcoa.scores;
   evals    = pcoa.evals;
   amongVar = evals./sum(evals);
   s         = size(scores,2); % # total eigenvalues
else % Constrained PCoA via db-RDA:
   %    rda      = f_rdaDB(dis,size(X,2),f_dummy(grps),0,0,0);
   %    scores   = rda.siteScores;
   %    amongVar = rda.evals./sum(rda.evals);
   %    s        = size(scores,2); % # non-zero canonical eigenvalues
   pcoa     = f_pcoa(dis,0,1,0); % scale eigenvectors, no negative eigenvalues
   scores   = pcoa.scores;
   evals    = pcoa.evals;
   scores   = scores(:,1:nClass-1);
   evals    = evals(1:nClass-1);
   s        = size(scores,2); % # non-zero canonical eigenvalues
   amongVar = evals./sum(evals);
end

% -----Create Canonical Plot:-----
if (nClass>1)
   figure('Name','RANDOM FOREST: Canonical Analysis');
else
   figure('Name','RANDOM FOREST: Unsupervised Learning');
end
set(gcf,'color','w');    % set bg color to white
hold on;

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
scores = (scores./repmat(max([max(abs(scores(:,1)));max(abs(scores(:,2)))]),...
   size(scores)));

hdl = zeros(1,noGrps); % preallocate
% plot points
for j = 1:noGrps
   gRows = find(grps==uGrps(j)); % get row indices for each group
   hdl(j) = plot(scores(gRows,1),scores(gRows,2),f_symb(j),'MarkerFaceColor',...
      f_rgb(j),'MarkerEdgeColor',f_rgb(j));
end

% Percent of total variation among GROUPS accounted for
% by Axis 1 & 2:
can1 = sprintf('%4.2f',(amongVar(1))*100);
can2 = sprintf('%4.2f',(amongVar(2))*100);

if (nClass>1)
   legend(hdl,Y_txt);
   
   xText = ['RF Canonical Axis I  (' num2str(can1) ' %)'];
   yText = ['RF Canonical Axis II (' num2str(can2) ' %)'];
   
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
else % Unsupervised learning:
   xText = ['RF Axis I  (' num2str(can1) ' %)'];
   yText = ['RF Axis II (' num2str(can2) ' %)'];
   
end
xlabel(xText,'FontWeight','bold');
ylabel(yText,'FontWeight','bold');

% Adjust figure properties:
axis tight;
axis(1.05*(axis)); % increase bounds of axis
axis equal
box on;
f_origin('hv'); % mark origin
f_arrow;        % re-draw arrows after changing axis scaling
% --------------------------------



% -----Biplot vectors:-----
ncX     = size(X,2);
corrX   = zeros(ncX,s); % preallocate
biplotX = zeros(ncX,s);

% Correlation with site scores:
for j=1:ncX
   for k=1:s
      corrX(j,k)   = f_corr(X(:,j),scores(:,k));     % correlation with site scores
      % biplotX(j,k) = corrX(j,k) * sqrt(amongVar(k)); % used for biplots
      biplotX = corrX; % do I want to use scaled or unscaled values?
   end
end

% Optionally trim biplot vectors:
biplotX = biplotX(idx,:);
X_txt   = X_txt(idx);

% Get length of each vector:
eDis = f_dis([0 0;biplotX(:,1:2)],'euc');
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
ylabel('Correlation with Canonical Axis II');

% Plot vectors:
for j = 1:size(biplotX,1)
   f_arrow([0 0],[biplotX(j,1) biplotX(j,2)],'size',[0.125]*0.75,'angle',20,...
      'Color','k');
   
   % Label X vectors:
   h = text(biplotX(j,1)*delta(j),biplotX(j,2)*delta(j),X_txt(j));
   set(h,'FontSize',8,'HorizontalAlignment','center','Color','k',...
      'Interpreter',fmt);
end

% Customize plot:
box on;
axis([-1 1 -1 1]*1.05)
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
f_origin('hv');
% --------------------------
