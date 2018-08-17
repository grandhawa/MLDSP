function f_RFplot(model,E,I,M)
% - create diagnostic plots of a Random Forest classifier
%
% USAGE: f_RFplot(model,E,I,M);
%
% E = plot out-of-box (OOB) error rates (= 1, default), include class rates (= 2)
% I = plot variable importance measures (default = 0)
% M = plot  margins of predictions (default = 0), unsorted (= 1), sorted (= 2)
% 
% SEE ALSO: f_RFclass, f_RFvis, f_RFregPlot

% -----Author:-----
% by David L. Jones, Sep-2009
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Oct-2011: updated documentation; check for regression vs. classification;

% -----Check input & set defaults:-----
if (nargin < 2), E = 1;    end % default plot OOB Error rates
if (nargin < 3), I = 0;    end % default don't plot importance measures
if (nargin < 4), M = 0;    end % default don't plot margins of predictions

if isequal(model.type,'regression')
   error('Use f_RFregPlot for REGRESSION!')
end

% Extract variables:
Y        = model.Y;
Y_txt    = model.Y_txt;
nClass   = model.nClass;
margin   = model.margin;
classAcc = model.classAcc;
meanAcc  = model.meanAcc;
gIndex   = model.gIndex;
errTr    = model.errtr;
N        = size(Y,1);

% Clean up:
clear model;

% Determine if variable importance measures were calculated:
if isnan(meanAcc(1))
   importance=0;
else
   importance=1;
end
   

% -----Create ERROR RATE plots:-----
if (E>0)
   figure('Name','OOB Error Rate');
   set(gcf,'color','w');    % set bg color to white
   box on;
   hold on;
   
   hdl = zeros(1,nClass+1); % preallocate
   
   % Plot average error rate:
   hdl(1) = plot(errTr(:,1),'color','k');
   
   if (E>1) % Plot class error rate:
      for j = 2:nClass+1
         hdl(j) = plot(errTr(:,j),'color',f_rgb(j-1),'LineStyle','--');
      end
      
      % Create legend:
      legend(hdl,['Mean'; Y_txt]);
   end
   
   % Set figure options:
   % title('OOB error rate');
   xlabel('\bf# Trees');
   yTxt = sprintf('Out of Box\nError Rate');
   ylabel(yTxt, 'FontWeight', 'bold');
   
   % Force Y-axis to start at 0:
   axisVar = axis;
   axis([axisVar(1) axisVar(2) 0 axisVar(4)]);
   
   % Use a gray-scale color map
   colormap(f_colormap(nClass));
end
% ----------------------------------




% -----PLOT IMPORTANCE MEASURES plots:-----
if (I>0) && (importance>0)
   figure('Name','Variable Importance');
   set(gcf,'color','w');    % set bg color to white
   
   subplot(2,1,1);
   bar(meanAcc,0.5);
   xlabel('\bfVariable');
   yTxt = sprintf('Mean decrease in\nAccuracy');
   ylabel(yTxt,'FontWeight','bold');
   title('\bfVariable Importance');
   
   subplot(2,1,2);
   bar(gIndex,0.5);
   xlabel('\bfVariable');
   yTxt = sprintf('Mean decrease in\nGini Index');
   ylabel(yTxt,'FontWeight','bold');
   
   figure('Name','Class Importance');
   set(gcf,'color','w');    % set bg color to white
   
   hdl = bar(classAcc','grouped'); % transpose so row=variable, col=class
   xlabel('\bfVariable');
   yTxt = sprintf('Mean decrease in\nAccuracy');
   ylabel(yTxt,'FontWeight','bold');
   titleTxt = sprintf('Class-specific Importance\n(legend refers to class)');
   title(titleTxt,'FontWeight','bold');
   
   % Create legend:
   legend(hdl,[Y_txt],'location','Best');
   
   % Use a gray-scale color map
   colormap(f_colormap(nClass));
end


% -----Plot MARGIN of PREDICTIONS:-----
if (M>0)
   figure('Name','Margin of Predictions');
   set(gcf,'color','w');    % set bg color to white
   box on;
   hold on;
   if M>1 % sort observations by margin:
      sMargin = sortrows([margin Y],1);
   else   % unsorted so X-axis indicates index to obs:
      sMargin = [margin Y]; 
   end
   for i = 1:N % plot classes with different colors
      plot(i,sMargin(i,1),'LineStyle','none','Marker','o','MarkerFaceColor',...
         f_rgb(sMargin(i,2)),'MarkerEdgeColor','none');
   end
   xlabel('\bfObservation');
   ylabel('\bfMagnitude');
   title('\bfMargin of Predictions');
end
% -------------------------------------
