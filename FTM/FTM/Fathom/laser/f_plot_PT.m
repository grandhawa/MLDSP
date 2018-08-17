function f_plot_PT(X,A,SRM,lod,log,sm,cps)
% - plot 'otolith profile' surface transect data (PPM)
%
% USAGE: f_plot_PT(X,A,SRM,LOD,log,sm,cps)
%
% X   = structure of 'otolith profile' surface transect data created by f_cps2ppm_PT
% A   = cell array specifying which analytes in X to plot; if empty all are potted
%       e.g., A = {'Li7' 'Sr88' 'U238'};
% SRM = name of variables in file SRM.mat to plot 95% confidence envelopes
%      e.g., SRM = {'nist610','nist612','nist614'};              (default = none)
% LOD = plot corresponding limits of detection                      (default = 0)
% log = use a log-scale for the Y-axis                              (default = 1)
% sm  = filter/smooth time-series data                              (default = 0)
% cps = plot raw cps data (= 1), 'cps/total cps' (= 2),
%       or molar ratios (= 3)                                       (default = 0)
%
% SEE ALSO: f_cps2ppm_PT, f_export_PT

% -----Notes:-----
% When using the LOG scale: ppm's = 1 are plotted as 0's while ppm's < 1
% are ignored and appear as holes in the timeseries. So, if values that were
% below the LOD's were set = 0, these values are NOT plotted in the time
% series.
% 
% If 'SM=1' the PPM data along with the CPS or MOLAR RATIO data are smoothed
% too.

% -----Author:-----
% by David L. Jones, July-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2011: added support for 95% confidence envelopes
% Oct-2011: updated documentation, added support for LOD's; log scale is
%           now optional; added support for linear scale
% Jan-2012: added support of filter/smooth data
% Mar-2012: added title and changed line style for LOD plot
% Jun-2012: adjust Y-axis label according to value of cps
% Jun-2013: now optionally plots molar ratios; now also smooths the CPS or Molar
%           Ratio data if SM=1

% TO DO: check for existence of *.pos field before continuing?

% -----Check input & set defaults:-----
if (nargin < 2), A     = []; end % default set empty to plot all
if (nargin < 3), SRM   = []; end % default don't plot confidence envelopes
if (nargin < 4), lod   = 0;  end % default don't plot LOD's
if (nargin < 5), log   = 1;  end % default plot on a log scale
if (nargin < 6), sm    = 0;  end % default don't filter/smooth times-series data
if (nargin < 7), cps   = 0;  end % default don't plot raw cps or molar ratios

% Make sure A & SRM are compatible:
if ~isempty(SRM) && numel(A)>1
   error('To plot confidence envelopes, A must specify only 1 analyte!');
end

% Specify log or linear scale:
if (log>0)
   scale = 'log';
else
   scale = 'linear';
end

% Check whether data have been filtered/smoothed:
if (sm>0) && isequal(X.sm,{'data filtered'})
   error('Check SM: data have already been filtered!')
end

% Check if raw CPS data are available:
if (cps>0) 
   if (isfield(X,'cps'==0))
      error('X does not have a CPS field!')
   end
   if (isnan(X.cps))
      error('X.cps contains no data!')
   end
   
   switch cps
      case 0 % do nothing
      case 1 % raw cps
         RAW = X.cps; % extract 'cps' from input
      case 2 % ratio of analyte cps:total cps
         RAW = X.cps ./ repmat(sum(X.cps,2),1,size(X.cps,2));
      case 3 % molar ratios
         RAW = X.ratio; % extract 'ratio' from input
      otherwise
         error('CPS must be 0, 1, or 2!');
   end
end
% ----------------------

% Extract fields from input:
txt_iso = X.txt_iso(:);
pos     = X.pos;
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

% Retain only selected analytes:
txt_iso = txt_iso(LOC);
ppm     = ppm(:,LOC);
LOD     = LOD(:,LOC);
if (cps>0)
   RAW = RAW(:,LOC);
end

% Optionally filter/smooth times-series data:
if (sm>0)
   ppm   = f_filterSinclair(ppm);
   LOD   = f_filterSinclair(LOD);
   if (cps>0)
      RAW = f_filterSinclair(RAW);
   end
end

% -----Extract ppm values for confidence envelopes:-----
if (~isempty(SRM))
   
   % Check that file 'SRM.mat' is in your path:
   if (exist('SRM.mat','file')~=2)
      error('''SRM.mat'' is not in Matlab''s search path!')
   end
   
   nCE = numel(SRM); % get # of confidence envelopes to plot
   
   % Extract ppm values of target analyte from each SRM:
   for i = 1:nCE
      
      % Check that the 'SRM' specified is present:
      if (size(whos('-file','SRM.mat',SRM{i}),1)<1)
         eval(['error(''The variable ' SRM{i} ' is not present within SRM.mat!'')'])
      end
      
      % Load SRM:
      eval(['load(''SRM.mat'', ''' SRM{i} ''')']); % load specified variable
      eval(['f_rename ' SRM{i} ' SRM_i'])          % rename loaded variable to 'SRM_i'
      
      % Trim/sort SRM to match target analyte
      a          = regexprep(A,'\d',''); % remove digits from A
      [null,loc] = ismember(a,SRM_i.txt);
      
      % Collect values:
      SRM_txt(i) = SRM_i.txt(loc(loc>0));
      SRM_ppm(i) = SRM_i.ppm(loc(loc>0));
      SRM_sd(i)  = SRM_i.sd(loc(loc>0));
      clear null loc;
   end
end
% -------------------------------------------------

% Create plot:
figure('Name',inputname(1));
set(gcf,'color','w'); % set background color to white
hold on;

% Plot each analyte separately:
nc      = size(ppm,2); % # analytes to plot
hdl(nc) = NaN;         % preallocate

for i = 1:nc
   if nc>10
      hdl(i) = plot(pos,ppm(:,i),'Color',f_rgb(i),'LineStyle',f_style(i));
   else
      hdl(i) = plot(pos,ppm(:,i),'Color',f_rgb(i),'LineStyle','-');
   end
end

% -----Optionally plot confidence envelopes:-----
if (~isempty(SRM))
   
   % Mean values:
   ave      = []; % initialize
   ave_SD_p = [];
   ave_SD_m = [];
   for i = 1:nCE
      % Mean:
      ave = [ ave; [min(pos) SRM_ppm(i)] ];
      ave = [ ave; [max(pos) SRM_ppm(i)] ];
      ave = [ ave; [NaN NaN]             ];
      
      % Plus 2 SD:
      ave_SD_p = [ ave_SD_p; [min(pos) SRM_ppm(i)+(2*SRM_sd(i))] ];
      ave_SD_p = [ ave_SD_p; [max(pos) SRM_ppm(i)+(2*SRM_sd(i))] ];
      ave_SD_p = [ ave_SD_p; [NaN NaN]             ];
      
      % Minus 2 SD:
      ave_SD_m = [ ave_SD_m; [min(pos) SRM_ppm(i)-(2*SRM_sd(i))] ];
      ave_SD_m = [ ave_SD_m; [max(pos) SRM_ppm(i)-(2*SRM_sd(i))] ];
      ave_SD_m = [ ave_SD_m; [NaN NaN]             ];
   end
   
   % Plot 95% confidence envelopes:
   plot(ave(:,1),ave(:,2),'Color',[1 1 1]*0.70,'LineStyle','-'); % mean
   plot(ave_SD_p(:,1),ave_SD_p(:,2),'Color',[1 1 1]*0.70,'LineStyle','--'); % plus 2 SD
   plot(ave_SD_m(:,1),ave_SD_m(:,2),'Color',[1 1 1]*0.70,'LineStyle','--'); % minus 2 SD
end
% -----------------------------------------------

% Create legend:
legend(hdl,txt_iso);

% Set up Y-axis label:
switch cps
   case 0
      yTxt = 'Concentration (ppm)';
   case 1
      yTxt = 'raw CPS';
   case 2
      yTxt = 'Ratio of CPS/Total CPS';
   case 3
      yTxt = 'Molar Ratio';
end

% Customize plot:
set(gca,'Yscale',scale,'YLimMode','manual'); % log or linear scale
xlabel('Position (um)')
ylabel(yTxt);
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
   nc      = size(LOD,2); % # analytes to plot
   hdl(nc) = NaN;         % preallocate
   
   for i = 1:nc
      if nc>10
         hdl(i) = plot(pos,LOD(:,i),'Color',f_rgb(i),'LineStyle',f_style(i));
      else
         hdl(i) = plot(pos,LOD(:,i),'Color',f_rgb(i),'LineStyle','-');
      end
   end
   
   % Create legend:
   legend(hdl,txt_iso);
   
   % Customize plot:
   set(gca,'Yscale',scale,'YLimMode','manual'); % log or linear scale
   xlabel('Position (um)')
   ylabel('Concentration (ppm)');
   title('Limits of Detection');
   box on;
   
   % Customize axes:
   yMinMax = get(gca,'YLim');
   % set(gca,'Ylim',[1 yMinMax(2)+0.1*yMinMax(2)]); % truncate at 1
   set(gca,'Ylim',[0 yMinMax(2)+0.1*yMinMax(2)]);   % go to 0
   
   % Create axes labels with thousands separators:
   set(gca,'YtickLabel',regexprep( cellstr(num2str(get(gca,'YTick')')),'\d{1,3}(?=(\d{3})+(?!\d))', '$&,'))
end
% --------------------------------

% -----Optionally plot raw CPS:-----
if (cps>0)
   % Create plot:
   if (cps==1)
      figure('Name',[inputname(1) ' Raw counts-per-second']);
   elseif (cps==2)
      figure('Name',[inputname(1) ' Ratio of CPS to total CPS']);
   elseif (cps==3)
      figure('Name',[inputname(1) ' Molar Ratio']);
   end
   set(gcf,'color','w'); % set background color to white
   hold on;
   
   % Plot each analyte separately:
   nc      = size(RAW,2); % # analytes to plot
   hdl(nc) = NaN;             % preallocate
   
   for i = 1:nc
      if nc>10
         hdl(i) = plot(pos,RAW(:,i),'Color',f_rgb(i),'LineStyle',f_style(i));
      else
         hdl(i) = plot(pos,RAW(:,i),'Color',f_rgb(i),'LineStyle','-');
      end
   end
   
   % Create legend:
   legend(hdl,txt_iso);
   
   % Customize plot:
   set(gca,'Yscale',scale,'YLimMode','manual'); % log or linear scale
   xlabel('Position (um)')
   if (cps==1)
      ylabel('Counts-per-second (CPS)');
      title('Raw CPS');
   elseif (cps==2)
      ylabel('Ratio analyte cps/total cps');
      title('CPS ratio');
   elseif (cps==3)
      ylabel('Molar Ratio');
      title('Molar Ratio');
   end
   box on;
   
   if (cps<2)
      % Create axes labels with thousands separators:
      set(gca,'YtickLabel',regexprep( cellstr(num2str(get(gca,'YTick')')),'\d{1,3}(?=(\d{3})+(?!\d))', '$&,'))
   end
end
% --------------------------------

