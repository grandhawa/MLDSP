function hdl = f_cpsPlot(S,bc,tot,fmt,idx)
% - plot LA-ICP-MS transient signal data
%
% USAGE: hdl = f_cpsPlot(S,bc,tot,fmt,idx);
%
% S   = structure of time series data created by f_importXL or f_cpsParse
% bc  = plot background corrected signals                          (default = 0)
% tot = plot total counts vs. individual analytes                  (default = 0)
% fmt = use long format for Y-axis labels                          (default = 1)
% idx = index to analtyes to plot                           (default = plot all)
%
% hdl = handle to axis of plot for customizing ledgend
%       e.g., legend(hdl,cellArray);
%
% SEE ALSO: f_importXL, f_cpsParse

% -----Author:-----
% by David L. Jones, Nov-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 2), bc  = 0;               end % default don't plot background corrected signals
if (nargin < 3), tot = 0;               end % default don't plot total counts
if (nargin < 4), fmt = 1;               end % default use long format for Y axis labels
if (nargin < 5), idx = 1:size(S.cps,2); end % default plot all analytes

% For background correction, make sure inputs include a BG field:
if ( (bc>0) &&  (~ismember('bg',fieldnames(S))) )
   error('S does not contain a background region!')
end

% If '.txt_iso' field is not present, create one:
if ~ismember('txt_iso',fieldnames(S))
   S.txt_iso = cellstr(strcat(S.txt',regexprep( cellstr(num2str(S.iso')),'\s', '') ))';
end
% ---------------------------------------

% Optional background correction:
if (bc>0)
   Na             = size(S.cps,1); % # sweeps in Unknown signal
   S.cps          = S.cps - repmat(mean(S.bg.cps),Na,1);
   S.cps(S.cps<0) = 0;             % set negative background-corrected count rates to 0
end

% Create plot;
figure('Name','LA-ICP-MS Time Series Plot');
set(gcf,'color','w'); % set bg color to white

if (tot>0)
   % Plot total counts:
   plot(S.s,sum(S.cps,2));
else
   % Plot each analyte separately:
   nc      = numel(idx); % # analytes to plot
   hdl(nc) = NaN;        % preallocate
   
   hold on;
   for i = 1:nc
      hdl(i) = plot(S.s,S.cps(:,idx(i)),'Color',f_rgb(i),'LineStyle',f_style(i));
   end
   hold off;
   
   % Create legend:
   legend(hdl,S.txt_iso{idx});
end

% Customize plot:
set(gca,'Yscale','log','YLimMode','auto');
xlabel('Seconds')
ylabel('Signal Intensity (cps)');
box on;

if (fmt>0)
   % Create axes labels with thousands separators:
   set(gca,'YtickLabel',regexprep( cellstr(num2str(get(gca,'YTick')')),'\d{1,3}(?=(\d{3})+(?!\d))', '$&,'))
end
