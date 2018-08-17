function [LK,ML,reef] = m_binnedMap(plt)
% - utility file for exampleBubbleMap to aggregate sampling sites into bins
% 
% USAGE: [LK,ML,reef] = m_binnedMap(plt);
% 
% plt = boolean (1 = create plot, 0 = no plot)

% by David L. Jones, Sep-2008

% Load the data:
load mpa.mat bnp;        % BNP park boundaries
load mangrove.mat grove; % MANGROVE
load reef.mat bnpSites;  % REEF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                LEEWARD KEY:                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select a value for the number of bins that more/less evenly distributes data:
nBins = 4; [n,c] = hist(grove.lat(grove.str==1),nBins);
% > [n' c']
% ans =
%           102       25.415
%            97       25.448
%            53       25.482
%            72       25.515
%
% Get distance between bin centers:
inc = c(2) - c(1);
%
% Get boundaries of each bin:
bnd = [min(grove.lat(grove.str==1)):inc:max(grove.lat(grove.str==1))]';
%
% Specify bin membership:
bin = zeros(size(grove.lat,1),1);          % initialize
bin(grove.lat(grove.str==1) <  bnd(2))= 1; % bin 1
bin(grove.lat(grove.str==1) >= bnd(2))= 2; % bin 2
bin(grove.lat(grove.str==1) >= bnd(3))= 3; % bin 3 + 4
%
% Check against output from 'hist' above:
% > [sum(bin==1) sum(bin==2) sum(bin==3)]'
% ans =
%    102
%     97
%    125
%
% Get the LON of a site whose LAT is closest to the bin's center:
lon = zeros(3,1); % initialize
for i=1:3
   [null,idx] = min(abs(c(i) - grove.lat(grove.str==1 & bin==i)));
   lonVar     = grove.lon(grove.str==1 & bin==i);
   lon(i)     = lonVar(idx(1));
end
%
% Wrap results up into a structure:
LK.lon = lon;
LK.lat = [c(1:2) mean(c(3:4))]';
LK.bin = bin;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 MAINLAND:                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select a value for the number of bins that more/less evenly distributes data:
nBins = 6; [n,c] = hist(grove.lat(grove.str==2),nBins);
% > [n' c']
% ans =
%           118        25.42
%           133       25.466
%           102       25.512
%            97       25.559
%            90       25.605
%           117       25.651
%
% Get distance between bin centers
inc = c(2) - c(1);
%
% Get boundaries of each bin:
bnd = [min(grove.lat(grove.str==2)):inc:max(grove.lat(grove.str==2))]';
%
% Specify bin membership:
bin = zeros(size(grove.lat,1),1);          % initialize
bin(grove.lat(grove.str==2) <  bnd(2))= 1; % bin 1
bin(grove.lat(grove.str==2) >= bnd(2))= 2; % bin 2
bin(grove.lat(grove.str==2) >= bnd(3))= 3; % bin 3
bin(grove.lat(grove.str==2) >= bnd(4))= 4; % bin 4
bin(grove.lat(grove.str==2) >= bnd(5))= 5; % bin 5
bin(grove.lat(grove.str==2) >= bnd(6))= 6; % bin 6
%
% Check against output from 'hist' above:
% > [sum(bin==1) sum(bin==2) sum(bin==3) sum(bin==4) sum(bin==5) sum(bin==6)]'
% ans =
%    118
%    133
%    102
%     97
%     90
%    117
%
% Get the LON of a site whose LAT is closest to the bin's center:
lon = zeros(6,1); % initialize
for i=1:6
   [null,idx] = min(abs(c(i) - grove.lat(grove.str==2 & bin==i)));
   lonVar     = grove.lon(grove.str==2 & bin==i);
   lon(i)     = lonVar(idx(1));
end
%
% Wrap results up into a structure:
ML.lon = lon(:);
ML.lat = c';
ML.bin = bin;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 REEF:                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select a value for the number of bins that more/less evenly distributes data:
nBins = 6; [n,c] = hist(bnpSites.lat,nBins);
% > [n' c']
% ans =
%           119       25.323
%           171       25.381
%           187       25.439
%           187       25.498
%           158       25.556
%           109       25.615
%
% Get distance between bin centers
inc = c(2) - c(1);
%
% Get boundaries of each bin:
bnd = [min(bnpSites.lat):inc:max(bnpSites.lat)]';
%
% Specify bin membership:
bin = zeros(size(bnpSites.lat,1),1); % initialize
bin(bnpSites.lat <  bnd(2))= 1;      % bin 1
bin(bnpSites.lat >= bnd(2))= 2;      % bin 2
bin(bnpSites.lat >= bnd(3))= 3;      % bin 3
bin(bnpSites.lat >= bnd(4))= 4;      % bin 4
bin(bnpSites.lat >= bnd(5))= 5;      % bin 5
bin(bnpSites.lat >= bnd(6))= 6;      % bin 6
%
% Check against output from 'hist' above:
% > [sum(bin==1) sum(bin==2) sum(bin==3) sum(bin==4) sum(bin==5) sum(bin==6)]'
% ans =
%    119
%    171
%    187
%    187
%    158
%    109
%
% Get average LON of sites withing each bin:
lon = zeros(6,1); % initialize
for i=1:6
   lon(i) = mean(bnpSites.lon(bin==i));
end
% 
% Wrap results up into a structure:
reef.lon = lon(:);
reef.lat = c';
reef.bin = bin;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  M_MAP:                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plt>0)

   % -----M_Map setup:-----
   figure;               % open new figure window
   set(gcf,'color','w'); % set bg color so printed/exported lakes will be white

   % Create a base map:
   m_proj('mercator','longitudes',[-80.375 -80.0625],'latitudes',[25.25 25.6875]);
   m_usercoast('fmri.mat','patch',[0.5 0.80 0.5],'edgecolor','none');

   % Biscayne National Park:
   m_line(bnp(:,1),bnp(:,2),'LineWidth',0.5,'LineStyle','-','Color',[1 1 1]*0.85);


   % Leeward Key:
   m_line(LK.lon,LK.lat,'marker','o','markersize',4, 'MarkerFaceColor','b',...
      'MarkerEdgeColor','none','LineStyle','none');

   % Mainland:
   m_line(ML.lon,ML.lat,'marker','o','markersize',4, 'MarkerFaceColor','r',...
      'MarkerEdgeColor','none','LineStyle','none');

   % Reef sites:
   m_line(reef.lon,reef.lat,'marker','o','markersize',4,...
      'MarkerFaceColor','k','MarkerEdgeColor','none','LineStyle','none');

   % -----Complete Base Map:-----
   m_grid('box','fancy','fontsize',8,'linestyle','none',...
      'xtick',(-80.375:[7.5/60]:-80.0625),'ytick',(25.25:[7.5/60]:25.6875));

end
