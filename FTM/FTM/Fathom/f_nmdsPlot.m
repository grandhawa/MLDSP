function [hdl,vec] = f_nmdsPlot(result,sLabels,Y,yLabels,offset,fmt,wt)
% - plot results from a nonmetric Multidimensional Scaling (nMDS)
%
% USAGE: [hdl,vec] = f_nmdsPlot(result,sLabels,Y,yLabels,offset,fmt,wt);
%
% result  = structure of results obtained from f_nmds
%
% sLabels = cell array of site labels; if empty, plot as filled circles
%           e.g., sLabels = {'A' 'A' 'A' 'B' 'B' 'C'};
%
% Y        = matrix of (transformed) response variables
% yLabels = cell array of Y labels; if empty, autocreate
%           e.g., yLabels = {'sp1' 'sp2' 'sp3'};
%
% offset  = yLabel offset factor ranging from 0-1                  (default = 0)
% fmt     = format of labels; 'tex' = TeX formatting (default) or 'none'
% wt      = weighted species biplot vectors                        (default = 0)
%
% hdl      = handle to axis of group plot for customizing ledgend
%            e.g., legend(hdl,cellArray);
% vec      = biplot vector coordinates
%
% SEE ALSO: f_nmds, f_pcoa

% -----Notes:-----
% When the response data represent species abundance data set 'wt=1' to
% generate weighted species biplot vectors, otherwise use 'wt=0' to
% generate standard correlation vectors.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp.
% Legendre, P. and E. Gallagher. Ecologically meaningful transformations
%   for ordination biplots of species data. Oecology 129: 271-280.

% -----Author:-----
% by David L. Jones, Mar-2013
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 2), sLabels = [];    end % default don't plot site labels
if (nargin < 3), Y       = [];    end % no biplot
if (nargin < 4), yLabels = [];    end % default Y labels
if (nargin < 5), offset  = 0;     end % set default offset to 0
if (nargin < 6), fmt     = 'tex'; end % default use TeX formatting
if (nargin < 7), wt      = 0;     end % default no weighting for biplot vectors

% Check offset:
if (offset<0) || (offset>1)
   error('Offset must be in the range 0-1!')
end

if ~isempty(Y)
   % Create/check yLabels:
   ncY = size(Y,2); % get # response variables:
   
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
end
% ---------------------------------------

% Unwrap variables from structure:
scores = result.scores;

% If only 1 nMDS axis, add 2nd jittered axis:
if (size(scores,2)==1)
   a      = min(scores) * 0.25;
   b      = max(scores) * 0.25;
   jit    = a + (b-a).*rand(size(scores));
   scores = [scores jit];
else
   jit = [];
end

% -----Create nMDS Plot:-----
figure('Name','nMDS: nonmetric Multidimensional Scaling');
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

title('\bfnonmetric Multidimensional Scaling');
xText = 'Dimension I';
if isempty(jit)
   yText = 'Dimension II';
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


% -----Biplot vectors:-----
if (~isempty(Y))
   s = size(scores,2);
   
   % Correlation of original Y variables with each nMDS axis:
   vec = zeros(ncY,s); % preallocate
   for i = 1:s
      for j = 1:ncY
         if (wt>0) % Weighted species biplot scores:
            vec(j,i) = f_corr(Y(:,j),scores(:,i))*std(Y(:,j))/std(scores(:,i));
         else
            vec(j,i) = f_corr(Y(:,j),scores(:,i));
         end
      end
   end
   
   if (wt>0)
      % Scale vectors to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
      vec = (vec./repmat(max([max(abs(vec(:,1)));max(abs(vec(:,2)))]), size(vec)));
   end
   
   % Get length of each vector:
   eDis = f_dis([0 0;vec(:,1:2)],'euc');
   eDis = eDis(2:end,1);
   
   % Set delta inversely proportional to length of longest vector:
   ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
   ratio = ones(size(ratio))./ratio;
   delta = 1 + (1* repmat(offset,size(eDis)) .* ratio);
   
   if (wt>0)
      figure('Name','Weighted Biplot Vectors');
   else
      figure('Name','Correlation Vectors');
   end
   set(gcf,'color','w'); % set bg color to white
   hold on;
   
   nCols  = 2;
   
   if (wt<1)
      % Draw circle for correlation vectors:
      radius = sqrt(2/nCols);
      rectangle('Position',[-radius,-radius,radius*2,radius*2],...
         'Curvature',[1,1],'EdgeColor',[1 1 1]*0.75);
   
      xText = 'Correlation with Dim I';
      yText = 'Correlation with Dim II';
   else
      xText = 'Dimension I';
      yText = 'Dimension II';
   end
   
   xlabel(xText);
   if isempty(jit)
      ylabel(yText);
   else
      ylabel('Jittered Axis');
   end
   
   axis([-1 1 -1 1]);
   daspect([1 1 1]);
   
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
   
else
   vec = NaN;
end
% ------------------------------------------
