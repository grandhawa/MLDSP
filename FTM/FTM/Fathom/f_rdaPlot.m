function H = f_rdaPlot(result,Y,wascores,scale,offset,yLabels,xLabels,fmt,iter,sLabels)
% - ordination distance biplot for a Redundancy Analysis (RDA)
%
% USAGE: H = f_rdaPlot(result,Y,wascores,scale,offset,yLabels,xLabels,'fmt',iter,sLabels)
%
% result   = structure of results obtained from f_rda or f_rdaDB
% Y        = matrix of (transformed) response variables
% wascores = plot Y's weighted average scores (vs. biplot vectors) (default = 1)
% scale    = scaling factor for biplot vectors                 (default = [1 1])
%            use [scaleX scaleY] for different scaling
% offset   = label offset                                          (default = 0)
%
% yLabels  = cell array of Y labels; if empty, autocreate
%            e.g., yLabels = {'sp1' 'sp2' 'sp3'};
%
% xLabels  = cell array of X labels; if empty, autocreate
%            e.g., xLabels = {'var1' 'var2' 'var3'};
% 
% fmt       = format of labels; 'tex' = TeX formatting (default) or 'none'
% iter      = iterations for permutation test of axis correlation  (default = 0)
% 
% sLabels  = cell array of site labels; if empty, plot as filled circles
%            e.g., sLabels = {'A' 'B' 'C'};
% 
% 
% H = structure with the following fields:
%  .biplotY    = optionally scaled species scores          (= vector endpoints)
%  .biplotX    = optionally scaled environmental scores    (= vector endpoints)
%  .p          = p-value of correlation of explanatory variables with frist 2
%                canonical axes 
%  .siteScores = site scores
%
% SEE ALSO: f_rda, f_rdaAIC, f_rdaStepwise, f_rdaDB

% -----Dependencies:-----
% This function requires A. J. Johnson's ARROW function, which is
% included in the Fathom Toolbox as f_arrow

% -----References:-----
% Legendre, P. and L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp.
% Legendre, P. and E. Gallagher. Ecologically meaningful transformations for
%   ordination biplots of species data. Oecology 129: 271-280. 
% McCune, B. 1997. Influence of noisy environmental data on canonical
%   correspondence analysis. Ecology 78: 2617-2623.

% -----Author:-----
% by David L. Jones, Oct-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: allow differential scaling of X,Y biplot vectors, fixed bug
%           in default label creation
% Aug-2003: set bg color to white, changed symbol/size for siteScores
% Feb-2008: updated docs; changed & to &&; biplotX is now created here vs.
%           imported from f_rda; changed num2Str to num2str, changed colors of
%           wascores and vectors
% Mar-2008: output Y,X biplot endpoints so, for example, you can re-create a
%           db-RDA plot with only those taxa having the strongest correlation
% Apr-2008: check size of labels; default labels are created for colums (vs.
%           rows) of input data; optionally format labels as TeX 
% Jun-2009: added ITER to allow significance testing of correlation of
%           explanatory variables with each canonical axis, output H
% Mar-2010: clarified documentation; call f_arrow after redrawing axes; center
%           response variable
% Nov-2011: don't center response variable (doesn't allow wascores)
% Jul-2012: now jitters Axis II if there's only one canonical axis (f_rdaDB
%           will produce canonical axes for univariate data, but not f_rda)
% Dec-2012: check input for FMT

% -----Unwrap variables from structure:-----
fitScores  = result.fitScores;
siteScores = result.siteScores;
X          = result.X;
canVar     = result.canVar;

% ----- Check input and set default values: -----
if (size(fitScores,1) == size(Y,1))<1
   error('The # of rows in Y & RESULT.FITSCORES must be equal!');
end;

if (nargin < 3), wascores = 1; end % plot Y's WA scores
if (nargin < 4), scale    = 1; end % set default scaling to 1
if (nargin < 5), offset   = 0; end % default is no label offset
if (nargin < 6), yLabels  = cellstr(num2str([1:size(Y,2)]'))'; end % default Y labels
if (nargin < 7), xLabels  = cellstr(num2str([1:size(X,2)]'))'; end % default X labels
if (nargin < 8), fmt      = 'tex'; end % default use TeX formatting
if (nargin < 9), iter     = 0;     end % default no significance test of correlations
if (nargin < 10), sLabels = [];    end % default no site labels

% If labels are not cell arrays, try forcing them:
if iscell(yLabels)<1, yLabels = num2cell(yLabels); end
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end
if ~isempty(sLabels) && iscell(sLabels)<1, sLabels = num2cell(sLabels); end

ncX = size(X,2); % # explanatory variables
ncY = size(Y,2); % # response variables

% Make sure labels are of compatible size:
xLabels = xLabels(:); % force cell array into a column
yLabels = yLabels(:);
if size(yLabels,1) ~= ncY
   error('Size of yLables doesn''t match # of Y variables!')
end
if size(xLabels,1) ~= ncX
   error('Size of xLables doesn''t match # of X variables!')
end

% Check FMT
if ( ~isequal(fmt,'tex') && ~isequal(fmt,'none') )
   error('FMT must be ''tex'' or ''none''!')
end
% -----------------------------------------------

% If only 1 canonical axis, add 2nd jittered axis:
if (size(siteScores,2)==1)
   a          = min(siteScores) * 0.25;
   b          = max(siteScores) * 0.25;
   c          = min(fitScores) * 0.25;
   d          = max(fitScores) * 0.25;
   jit        = a + (b-a).*rand(size(siteScores));
   jit_2      = c + (d-c).*rand(size(fitScores));
   siteScores = [siteScores jit];
   fitScores  = [fitScores jit_2];
   canVar     = [canVar;1];
else
   jit = [];
end

crds = fitScores(:,1:2);
s    = size(siteScores,2); % # non-zero canonical eigenvalues

% Response/Explanatory variables:
X = f_stnd(X);   % standardize, especially important if measured on different scales

% -----Correlations with original variables:-----
% The square of these correlations give the fraction of variance of each
% variable expressed along each canonical axis. 

corrX     = zeros(ncX,s); % preallocate
corrX_fit = zeros(ncX,s); 
biplotX   = zeros(ncX,s); 

for j=1:ncX
   for k=1:s
      corrX(j,k)     = f_corr(X(:,j),siteScores(:,k));    % correlation with site scores
      corrX_fit(j,k) = f_corr(X(:,j),fitScores(:,k));     % correlation with fitted site scores
      biplotX(j,k)   = corrX_fit(j,k) * sqrt(canVar(k));  % used for biplots
   end
end

% -----Tabulate p-values:-----
% Some of the explanatory variables may not be significantly correlated with the
% plotted canonical axes (Axis I & II), so we can calculate the p-values of
% these correlations for each of the explanatory variables:

if iter>0
   p         = zeros(ncX,2); % preallocate, only interested in first 2 axes
   counter   = 0;
   for j=1:ncX
      for k = 1:2
         counter = counter+1;
         fprintf('\nCalculating p-value for variable %d of %d...\n',counter,ncX*2);
         [null p(j,k)] = f_corr(X(:,j),fitScores(:,k),0,iter); % correlation with site scores
      end
   end
else
   p = NaN;
end


% -----This was cut from F_RDA:-----
% % If Y = PCoA axes derived from a distance matrix (e.g., from f_rdaDB)
% % you would want to calculate these correlations from the original (transformed)
% % data used to create the distance matrix
% %
% corrY     = zeros(ncY,s); % preallocate
% corrY_fit = zeros(ncY,s);
% 
% for j=1:ncY
%    for k=1:s
%       corrY(j,k)     = f_corr(Y(:,j),siteScores(:,k));
%       corrY_fit(j,k) = f_corr(Y(:,j),fitScores(:,k));
%    end
% end
% ----------------------------------
% -----------------------------------------------


% =========================================================================
std_crds = std(crds);    % get standard deviations
std_Y    = std(Y); 

biplotY  = zeros(ncY,2); % preallocate

for i = 1:ncY
   biplotY(i,1) = [f_corr(Y(:,i),crds(:,1))]*std_Y(:,i)/std_crds(1);  % 1st axis
   biplotY(i,2) = [f_corr(Y(:,i),crds(:,2))]*std_Y(:,i)/std_crds(2);  % 2nd axis
end

% -----Scale biplots:-----
scale = scale(:)';     % force column vector

if (size(scale,2) < 2) % scale vectors in same way if only 1 value given
   scale = [scale scale];
end

scaleX = scale(1); % scale biplot vectors of explanatory variables
scaleY = scale(2); % scale biplot vectors of response variable

biplotX = biplotX*scaleX;
biplotY = biplotY*scaleY;

% =========================================================================

% -----Create Canonical Plot:-----
figure;
title('\bfRedundancy Analysis');
set(gcf,'color','w'); % set bg color to white
hold on;   

% Plot siteScores:
if isempty(sLabels)
   plot(siteScores(:,1),siteScores(:,2),'o','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
else
   plot(siteScores(:,1),siteScores(:,2),'.w','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
   text(siteScores(:,1),siteScores(:,2),sLabels,'HorizontalAlignment','center',...
      'VerticalAlignment','middle','Color','b');
end

if (wascores>0) % plot Y's WA scores (McCune, 1997):
   waY = f_wascores(siteScores,Y,0);
   h   = text(waY(:,1),waY(:,2),yLabels);
   set(h,'HorizontalAlignment','center','Color',f_rgb(2),'FontSize',8,...
      'FontWeight','bold', 'Interpreter', fmt);      
   
else % Plot Y biplot vectors:
   for j = 1:ncY
      f_arrow([0 0],[biplotY(j,1) biplotY(j,2)],'size',[0.125]*1,'angle',20,...
         'Color',f_rgb(2));
      
      if (biplotY(j,1)>=0 && biplotY(j,2)>=0)
         delta1 = offset;
         delta2 = offset;
      elseif (biplotY(j,1)<0 && biplotY(j,2)>=0)     
         delta1 = -1*offset;
         delta2 =    offset;
      elseif (biplotY(j,1)<0 && biplotY(j,2)<0)   
         delta1 = -1*offset;
         delta2 = -1*offset;
      elseif (biplotY(j,1)>=0 && biplotY(j,2)<0)
         delta1 =    offset;
         delta2 = -1*offset;
      end
      
      % Label Y vectors:
      h = text(biplotY(j,1)+(delta1*0.5),biplotY(j,2)+(delta2*0.5),yLabels(j));
      set(h,'FontSize',8,'HorizontalAlignment','center','Color',f_rgb(2),...
         'Interpreter', fmt);
   end
end

% -----Plot X vectors:-----
for j = 1:size(biplotX,1)
   f_arrow([0 0],[biplotX(j,1) biplotX(j,2)],'size',[0.125]*0.75,'angle',20,...
      'Color','r');
   
   if (biplotX(j,1)>=0 && biplotX(j,2)>=0)
      delta1 = offset;
      delta2 = offset;
   elseif (biplotX(j,1)<0 && biplotX(j,2)>=0)     
      delta1 = -1*offset;
      delta2 =    offset;
   elseif (biplotX(j,1)<0 && biplotX(j,2)<0)   
      delta1 = -1*offset;
      delta2 = -1*offset;
   elseif (biplotX(j,1)>=0 && biplotX(j,2)<0)
      delta1 =    offset;
      delta2 = -1*offset;
   end;
   
   % Label X vectors:
   h = text(biplotX(j,1)+(delta1*0.5),biplotX(j,2)+(delta2*0.5),xLabels(j));
   set(h,'FontSize',8,'HorizontalAlignment','center','Color','r',...
      'Interpreter',fmt);
end;
%--------------------------

% Percent of total variation in Y accounted for by Axis 1 & 2:
can1 = sprintf('%2.2f',canVar(1)*100);
can2 = sprintf('%2.2f',canVar(2)*100);
xText = ['Canonical Axis I (' num2str(can1) ' %)'];
if isempty(jit)
   yText = ['Canonical Axis II (' num2str(can2) ' %)'];   
else
   yText = 'Jittered Axis';
end
xlabel(xText);
ylabel(yText);   

% Set figure options:
axis tight;
axis(1.05*(axis));
axis equal;
box on;

% Redraw arrows after axes changes:
f_arrow;

% Mark the origin:
f_origin('hv',':');

%-----Wrap results up into a structure:-----
% You may need this in order to re-create the plot with only a subset of, say,
% the Y data
H.biplotY    = biplotY;
H.biplotX    = biplotX;
H.p          = p;
H.siteScores = siteScores;

