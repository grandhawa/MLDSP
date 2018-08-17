function f_plot_PS(X,A,lod,log)
% - plot 'otolith profile' data collected as a series of spots
%
% USAGE: f_plot_PS(X,A,LOD,log)
%
% X   = structure of 'otolith spot profile' data created by f_cps2ppm_PS
% A   = cell array specifying which analytes in X to plot; if empty all are potted
%       e.g., A = {'Li7' 'Sr88' 'U238'};
% LOD = plot corresponding limits of detection                      (default = 0)
% log = use a log-scale for the Y-axis                              (default = 1)
%
% SEE ALSO: f_cps2ppm_PS

% -----Notes:-----
% When using the LOG scale: ppm's = 1 are plotted as 0's while ppm's < 1
% are ignored and appear as holes in the timeseries. So, if values that were
% below the LOD's were set = 0, these values are NOT plotted in the time
% series.

% -----Author:-----
% by David L. Jones, Oct-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: updated documentation

% -----Check input & set defaults:-----
if (nargin < 2), A    = []; end % default set empty to plot all
if (nargin < 3), lod  = 0;  end % default don't plot LOD's
if (nargin < 4), log  = 1;  end % default plot on a log scale

% Specify log or linear scale:
if (log>0)
   scale = 'log';
else
   scale = 'linear';
end
% ----------------------

% Extract fields from input:
txt_iso = X.txt_iso(:);
spot    = X.spot;
ppm     = X.ppm;
LOD     = X.LOD;

% Check list of analytes:
if isempty(A)
   A = txt_iso; % plot all analytes if none are specified
end
if ~iscell(A), error('A should be a cell array of strings!'); end
A = A(:); % force column vector

% Get index of analytes to plot:
[TF,LOC] = ismember(A,txt_iso);
idx = find(TF==0);
if ~isempty(idx)
   fprintf('These analytes were not found in %s:\n',inputname(1));
   fprintf('  %s \n\n',A{idx})
   error('All analytes in A must be found in X!');
end

% Create plot:
figure('Name',inputname(1));
set(gcf,'color','w'); % set background color to white
hold on;

% Plot each analyte separately:
nc      = numel(LOC); % # analytes to plot
hdl(nc) = NaN;        % preallocate

for i = 1:nc
   if nc>10
      hdl(i) = plot(spot,ppm(:,LOC(i)),'Color',f_rgb(i),'LineStyle',f_style(i),...
         'Marker',f_symb(i));
   else
      hdl(i) = plot(spot,ppm(:,LOC(i)),'Color',f_rgb(i),'LineStyle','-','Marker','.');
   end
end

% Create legend:
legend(hdl,txt_iso{LOC});

% Customize plot:
set(gca,'Yscale',scale,'YLimMode','manual'); % log or linear scale
xlabel('Spot #')
ylabel('Concentration (ppm)');
box on;

% Customize axes:
yMinMax = get(gca,'YLim');
% set(gca,'Ylim',[1 yMinMax(2)+0.1*yMinMax(2)]); % truncate at 1
set(gca,'Ylim',[0 yMinMax(2)+0.1*yMinMax(2)]);   % go to 0

% Create axes labels with thousands separators:
set(gca,'YtickLabel',regexprep( cellstr(num2str(get(gca,'YTick')')),'\d{1,3}(?=(\d{3})+(?!\d))', '$&,'))


% -----Optionally plot LOD's:-----
if lod>0
   % Create plot:
   figure('Name',[inputname(1) ' Limits of detection']);
   set(gcf,'color','w'); % set background color to white
   hold on;
   
   % Plot each analyte separately:
   nc      = numel(LOC); % # analytes to plot
   hdl(nc) = NaN;        % preallocate
   
   for i = 1:nc
      if nc>10
         hdl(i) = plot(spot,LOD(:,LOC(i)),'Color',f_rgb(i),'LineStyle',f_style(i),...
            'Marker',f_symb(i));
      else
         hdl(i) = plot(spot,LOD(:,LOC(i)),'Color',f_rgb(i),'LineStyle',':', 'Marker','.');
      end
   end
   
   % Create legend:
   legend(hdl,txt_iso{LOC});
   
   % Customize plot:
   set(gca,'Yscale',scale,'YLimMode','manual'); % log or linear scale
   xlabel('Position (um)')
   ylabel('Concentration (ppm)');
   box on;
   
   % Customize axes:
   yMinMax = get(gca,'YLim');
   % set(gca,'Ylim',[1 yMinMax(2)+0.1*yMinMax(2)]); % truncate at 1
   set(gca,'Ylim',[0 yMinMax(2)+0.1*yMinMax(2)]);   % go to 0
   
   % Create axes labels with thousands separators:
   set(gca,'YtickLabel',regexprep( cellstr(num2str(get(gca,'YTick')')),'\d{1,3}(?=(\d{3})+(?!\d))', '$&,'))
end
% --------------------------------
