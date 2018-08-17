function hdl = f_pcaPlot(result,sLabels,yLabels,offset,fmt)
% - plot results from a Principal Components Analysis
%
% USAGE: hdl = f_pcaPlot(result,sLabels,yLabels,offset,fmt);
%
% result  = structure of results obtained from f_pca
%
% sLabels = cell array of site labels; if empty, plot as filled circles
%           e.g., sLabels = {'A' 'A' 'A' 'B' 'B' 'C'};
%
% yLabels = cell array of Y labels; if empty, autocreate
%           e.g., yLabels = {'env1' 'env2' 'env3'};
%
% offset  = yLabel offset factor ranging from 0-1                  (default = 0)
% fmt     = format of labels; 'tex' = TeX formatting (default) or 'none'
%
% hdl      = handle to axis of group plot for customizing ledgend
%            e.g., legend(hdl,cellArray);
%
% SEE ALSO: f_pca, f_pcoaPlot

% -----Notes:-----
% The EQUILIBRIUM CIRCLE provides a graphical method of determining the
% relative contribution each VARIABLE makes to the formation of the reduced
% spaced defined by PC Axis I & II. Only VARIABLES whose vectors extend to or
% beyond the radius of the Equilibrium Circle are considered to have made a
% significant contribution to that reduced space.

% -----References:-----
% Elmore, K. L., and M. B. Richman. 2001. Euclidean distance as a
%   similarity metric for principal component analysis. Monthly Weather
%   Review 129: 540-549.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, May-2012
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% Jan-2013: fixed problem with autocreating yLabels

% -----Set defaults and check input:-----
if (nargin < 2), sLabels = [];    end % default don't plot site labels
if (nargin < 3), yLabels = [];    end % default Y labels
if (nargin < 4), offset  = 0;     end % set default offset to 0
if (nargin < 5), fmt     = 'tex'; end % default use TeX formatting

% Check offset:
if (offset<0) || (offset>1)
   error('Offset must be in the range 0-1!')
end

% Unwrap variables from structure:
scores = result.scores;
vec    = result.evects;
expl   = result.expl;
txt    = result.txt;

% Autocreate yLabels:
if isempty(yLabels)
   yLabels = cellstr(num2str((1:size(vec,2))'))';
end

% If labels are not cell arrays, try forcing them:
if iscell(yLabels)<1, yLabels = num2cell(yLabels); end

% Make sure labels are of compatible size:
yLabels = yLabels(:);
if size(yLabels,1) ~= size(vec,2)
   error('Size mismatch b/n yLables and VEC!')
end
% ---------------------------------------

% If only 1 PCA axis, add 2nd jittered axis:
if (size(scores,2)==1)
   a      = min(scores) * 0.25;
   b      = max(scores) * 0.25;
   jit    = a + (b-a).*rand(size(scores));
   scores = [scores jit];
else
   jit = [];
end

% -----Create PCA Plot:-----
figure('Name','PCA: Principal Components Analysis');
set(gcf,'color','w');    % set bg color to white
hold on;

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
scores = (scores./repmat(max([max(abs(scores(:,1)));max(abs(scores(:,2)))]), size(scores)));

% Plot sites:
if isempty(sLabels)
   hdl = plot(scores(:,1),scores(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
else
   hdl = plot(scores(:,1),scores(:,2),'.w','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
   text(scores(:,1),scores(:,2),sLabels,'HorizontalAlignment','center',...
      'VerticalAlignment','middle','Color','b','Interpreter',fmt);
end

% Percent of total variation accounted for:
axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));

title(txt);
xText = ['Axis I (' num2str(axis1) ' %)'];
if isempty(jit)
   yText = ['Axis II (' num2str(axis2) ' %)'];
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


% -----Equilibrium Circle:-----
figure('Name','Equilibrium Circle');
set(gcf,'color','w'); % set bg color to white
box on;
hold on;
hold on;
nCols  = size(vec,2);
radius = sqrt(2/nCols);
rectangle('Position',[-radius,-radius,radius*2,radius*2],...
   'Curvature',[1,1],'EdgeColor',[1 1 1]*0.75);
axis([-1 1 -1 1]);
daspect([1 1 1]);

title('Equilibrium Circle');
xText = 'PC Axis I';
yText = 'PC Axis II';
xlabel(xText);
ylabel('PC Axis II');
if isempty(jit)
   ylabel(yText);
else
   ylabel('Jittered Axis');
end

% Get length of each vector:
eDis = f_dis([0 0;vec(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio = ones(size(ratio))./ratio;
delta = 1 + (1* repmat(offset,size(eDis)) .* ratio);

n = size(vec,1); % get # of variables
for j = 1:n
   % Plot vectors:
   f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
   
   % Label vectors:
   h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),yLabels(j));
   set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter',fmt);
end

% Customize plot:
box on;
axis([-1 1 -1 1]*1.05)
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
f_origin('hv');
% -----------------------------

% 
% % Equilibrium Circle:
%    if (plt>2)
%       figure;
%       set(gcf,'color','w'); % set bg color to white
%       box on;
%       hold on;
%       radius = sqrt(2/nCols);
%       rectangle('Position',[-radius,-radius,radius*2,radius*2],...
%          'Curvature',[1,1],'EdgeColor',[0 0 1]);
%       axis([-1.1 1.1 -1.1 1.1]);
%       daspect([1 1 1]);
%       title('Equilibrium Circle');
%       xlabel('PC Axis I');
%       ylabel('PC Axis II');
%       % plot the vectors:
%       f_biplot(evects,1,0);
%       f_origin('hv');
%    end;
%    hold off



%    % Principal Component Scores:
%    figure;                             % opens new figure window
%    set(gcf,'color','w');               % set bg color to white
%    plot(scores(:,1),scores(:,2),'bo'); % plots column 1 x 2 of PC Scores
%    
%    pc1 = sprintf('%2.2f',(evals(1)/sum(evals))*100);
%    pc2 = sprintf('%2.2f',(evals(2)/sum(evals))*100);
%    xText = ['PC Axis I (' num2str(pc1) ' %)'];
%    yText = ['PC Axis II (' num2str(pc2) ' %)'];
%    xlabel(xText);
%    ylabel(yText);
%    f_origin('hv');



