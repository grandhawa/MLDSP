function f_CCorAplot(result,Y,yLabels,X,xLabels,sLabels,offset,fmt)
% - plot results from a CCorA analysis
%
% USAGE: f_CCorAplot(result,Y,yLabels,X,xLabels,sLabels,offset,fmt);
%
% result  = structure of results obtained from f_CCorA
% Y       = matrix of original (transformed) Y variables
% yLabels = cell array of Y labels; if empty, autocreate
%           e.g., yLabels = {'sp1' 'sp2' 'sp3'};
% X       = matrix of original (transformed) X variables
% yLabels = cell array of Y labels; if empty, autocreate
%           e.g., yLabels = {'temp' 'sal' 'depth'};
% sLabels = cell array of site labels; e.g., sLabels = {'A' 'B' 'C'};
%           []: if empty, plot as filled circles
%            0: don't plot sites
% offset  = label offset factor ranging from 0-1                   (default = 0)
% fmt     = format of labels; 'tex' = TeX formatting (default) or 'none'
%
% SEE ALSO: f_CCorA

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%  Elsevier Science BV, Amsterdam. xv + 853 pp.

% -----Author:-----
% by David L. Jones, Mar-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 6), sLabels = [];   end % default no site labels
if (nargin < 7), offset  = 0;    end % set default offset to 0
if (nargin < 8), fmt    = 'tex'; end % default use TeX formatting

% Check offset:
if (offset<0) || (offset>1)
   error('Offset must be in the range 0-1!')
end
% ---------------------------------------

% -----Unwrap variables from structure:-----
yScores = result.yScores;
xScores = result.xScores;

% -----Create/check Y Labels:-----
ncY = size(Y,2); % # Y variables:

% Autocreate yLabels:
if isempty(yLabels)
   yLabels = cellstr(num2str((1:ncY)'))';
end

% If labels are not cell arrays, try forcing them:
if iscell(yLabels)<1, yLabels = num2cell(yLabels); end

% Make sure labels are of compatible size:
yLabels = yLabels(:);
if size(yLabels,1) ~= ncY
   error('Size mismatch of yLables and Y!')
end

% -----Create/check X Labels:-----
ncX = size(X,2); % # X variables:

% Autocreate yLabels:
if isempty(xLabels)
   xLabels = cellstr(num2str((1:ncX)'))';
end

% If labels are not cell arrays, try forcing them:
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end

% Make sure labels are of compatible size:
xLabels = xLabels(:);
if size(xLabels,1) ~= ncX
   error('Size mismatch of xLables and X!')
end

% -----Create/check Site Labels:-----
nS = size(yScores,1); % # sites

% Autocreate yLabels:
if isempty(sLabels)
   sLabels = cellstr(num2str((1:size(X,1))'))';
end

% If labels are not cell arrays, try forcing them:
if ~isequal(0,sLabels) && iscell(sLabels)<1, sLabels = num2cell(sLabels); end

% Make sure labels are of compatible size:
if ~isequal(0,sLabels)
   sLabels = sLabels(:);
   if size(sLabels,1) ~= nS
      error('Size mismatch b/n sLabels and yScores!')
   end
end

% -----Create plot:-----
figure('Name','CCorA');
set(gcf,'color','w'); % set bg color to white

% -----Plot Y correlation vectors:-----
subplot(1,2,1);
hold on;
daspect([1 1 1]); % call before f_arrow

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
yScores = (yScores./repmat(max([max(abs(yScores(:,1)));max(abs(yScores(:,2)))]),...
   size(yScores)));

% Optionally plot sites:
if isequal(0,sLabels)
   plot(yScores(:,1),yScores(:,2),'o','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
else
   plot(yScores(:,1),yScores(:,2),'.w','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
   text(yScores(:,1),yScores(:,2),sLabels,'HorizontalAlignment','center',...
      'VerticalAlignment','middle','Color','b');
end

% Correlation of original Y variables with each canonical axis:
s = size(yScores,2);
yCor = zeros(ncY,s); % preallocate
for i = 1:s
   for j = 1:ncY
      yCor(j,i) = f_corr(Y(:,j),yScores(:,i));
   end
end

% Get length of each vector:
eDis = f_dis([0 0;yCor(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio = ones(size(ratio))./ratio;
delta = 1 + (1* repmat(offset,size(eDis)) .* ratio);

% Plot Y correlation vectors:
nCols = 2;
radius = sqrt(2/nCols);
rectangle('Position',[-radius,-radius,radius*2,radius*2],...
   'Curvature',[1,1],'EdgeColor',[1 1 1]*0.75);

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
title('\bfCCorA biplot (Y)');
xlabel('Canonical Axis I');
ylabel('Canonical Axis II');
axis([-1 1 -1 1]*1.05)
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
f_origin('hv',':',0.75);
f_arrow; % redraw arrows to avoid skewing
% ------------------------------------------


% -----Plot X correlation vectors:-----
subplot(1,2,2);
hold on;
daspect([1 1 1]); % call before f_arrow

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
xScores = (xScores./repmat(max([max(abs(xScores(:,1)));max(abs(xScores(:,2)))]),...
   size(xScores)));

% Optionally plot sites:
if isequal(0,sLabels)
   plot(xScores(:,1),xScores(:,2),'o','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
else
   plot(xScores(:,1),xScores(:,2),'.w','MarkerFaceColor',0.70*[1 1 1],...
      'MarkerEdgeColor','none');
   text(xScores(:,1),xScores(:,2),sLabels,'HorizontalAlignment','center',...
      'VerticalAlignment','middle','Color','r');
end

% Correlation of original X variables with each canonical axis:
s = size(xScores,2);
xCor = zeros(ncX,s); % preallocate
for i = 1:s
   for j = 1:ncX
      xCor(j,i) = f_corr(X(:,j),xScores(:,i));
   end
end

% Get length of each vector:
eDis = f_dis([0 0;xCor(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio  = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio  = ones(size(ratio))./ratio;
delta  = 1 + (1* repmat(offset,size(eDis)) .* ratio);

% Plot X correlation vectors:
nCols = 2;
radius = sqrt(2/nCols);
rectangle('Position',[-radius,-radius,radius*2,radius*2],...
   'Curvature',[1,1],'EdgeColor',[1 1 1]*0.75);
n = size(xCor,1); % get # of variables
for j = 1:n
   % Plot vectors:
   f_arrow([0 0],[xCor(j,1) xCor(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
   
   % Label vectors:
   h = text(xCor(j,1)*delta(j),xCor(j,2)*delta(j),xLabels(j));
   set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter',fmt);
end

% Customize plot:
box on;
title('\bfCCorA biplot (X)');
xlabel('Canonical Axis I');
ylabel('');
axis([-1 1 -1 1]*1.05)
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
f_origin('hv',':',0.75);
f_arrow; % redraw arrows to avoid skewing
% ------------------------------------------
