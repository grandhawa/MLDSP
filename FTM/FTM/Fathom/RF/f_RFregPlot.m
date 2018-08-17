function f_RFregPlot(model,E,I,norm)
% - create diagnostic plots of a Random Forest regression
%
% USAGE: f_RFregPlot(model,E,I,norm);
%
% E     = plot out-of-box (OOB) error rates  (default = 1)
% I     = plot variable importance measures  (default = 0)
% norm  = normalize importance measure by SD (default = 1)
% 
% SEE ALSO: f_RFreg, f_RFvis

% -----Notes:-----
% modified after f_RFclassPlot

% -----Author:-----
% by David L. Jones, Feb-2010
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Oct-2011: updated documentation; check for regression vs. classification;

% -----Check input & set defaults:-----
if (nargin < 2), E    = 1; end % default plot OOB Error rates
if (nargin < 3), I    = 0; end % default don't plot importance measures
if (nargin < 4), norm = 1; end % default don't normalize importance measures

if isequal(model.type,'classification')
   error('Use f_RFplot for CLASSIFICATION!')
end

% Extract variables:
meanAcc = model.meanAcc; % meanAcc (raw)
zScore  = model.zScore;  % meanAcc normalized by SD
meanMSE = model.meanMSE;
mse     = model.mse;

% Clean up:
clear model;

% Determine if variable importance measures were calculated:
if isnan(zScore(1))
   importance=0;
else
   importance=1;
end
   

% -----Create ERROR RATE plots:-----
if (E>0)
   figure('Name','OOB MSE');
   set(gcf,'color','w');    % set bg color to white
   box on;
   hold on;
         
   % Plot average error rate:
   plot(mse(:,1),'color','k');
        
   % Set figure options:
   xTxt = sprintf('# Trees');
   xlabel(xTxt, 'FontWeight', 'bold');
   
   yTxt = sprintf('Out-of-Box Error Rate');
   ylabel(yTxt, 'FontWeight', 'bold');
   
   % Force Y-axis to start at 0:
   axisVar = axis;
   axis([axisVar(1) axisVar(2) 0 axisVar(4)]);
end
% ----------------------------------

% -----PLOT IMPORTANCE MEASURES plots:-----
if (I>0) && (importance>0)
   figure('Name','Variable Importance');
   set(gcf,'color','w');    % set bg color to white
   
   subplot(2,1,1);
   if (norm>0)
      bar(zScore,0.5);
      yTxt = sprintf('%% Increase in MSE\n(normalized)');
   else
      bar(meanAcc,0.5);
      yTxt = sprintf('%% Increase in MSE\n(raw)');
   end
   xTxt = sprintf('Variable');
   xlabel(xTxt,'FontWeight','bold');
   ylabel(yTxt,'FontWeight','bold');
   
   subplot(2,1,2);
   bar(meanMSE,0.5);
   xlabel('\bfVariable');
   yTxt = sprintf('Decrease in \nNode Impurity (RSS)');
   ylabel(yTxt,'FontWeight','bold');
end
