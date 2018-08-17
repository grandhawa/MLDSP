% Example of creating bubble maps for a rather complex data set
% 
% by David L. Jones, Feb-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Load the data:
load mpa.mat bnp;            % BNP park boundaries
load mangrove.mat bay grove; % MANGROVE
load reef.mat bayReef;       % REEF

% Aggregate sites into bins:
[LK,ML,reef] = m_binnedMap(0);

% -----Setup annotations:-----
name_txt = {'Sergeant major' 'Yellowfin mojarra' 'French grunt'...
   'Sailor''s choice' 'Bluestriped grunt' 'Pinfish' 'Schoolmaster'...
   'Gray snapper' 'Rainbow parrotfish' 'Great barracuda'}';

stage_txt = {'juvenile' 'subadult' 'adult'}';
% ----------------------------

% Set plot options:
fnt = 6;  % set font size
mx  = 20; % set max symbol size

% Select which taxon to plot:
k = 1;

% -----M_Map setup:-----
figure;
set(gcf,'color','w'); % set bg color so printed/exported lakes will be white

maxM = zeros(3,1); % initialize max Mangrove
maxR = zeros(3,1); % initialize max Reef

% Initialize column pointer (there are 3 stages of 10 species)
col = k*3-2;

for i=1:3

   % LEEWARD KEY (average within bins, across years):
   meanLK{i} = zeros(3,1);
   for j=1:3
      idx          = find( (grove.str==1 & LK.bin==j) == 1);
      meanLK{i}(j) = mean(bay.numS(idx,col));
   end

   % MAINLAND (average within bins, across years):
   meanML{i} = zeros(6,1);
   for j=1:6
      idx          = find( (grove.str==2 & LK.bin==j) == 1);
      meanML{i}(j) = mean(bay.numS(idx,col));
   end

   % REEF (average within bins, across years):
   meanReef{i} = zeros(6,1);
   for j=1:6
      idx         = find( (reef.bin==j) == 1);
      meanReef{i}(j) = mean(bayReef.numS(idx,col));
   end

   maxM(i) = max([meanLK{i};meanML{i}]);
   maxR(i) = max(meanReef{i});

   col = col+1; % increment to next stage
end

% Scale the maximum symbol in each subplot according to the maximum abundance
% across plots (separately for Mangrove and Reef):
scaleM = maxM / max(maxM);
scaleR = maxR / max(maxR);


for i = 1:3 % Create a Bubble plot for each of 3 stages:
   subplot(1,3,i);

   % Create a base map:
   m_proj('mercator','longitudes',[-80.375 -80.0625],'latitudes',[25.25 25.6875]);
   m_usercoast('fmri.mat','patch',[0.5 0.80 0.5],'edgecolor','none');
   % m_usercoast('bnp_h.mat','patch',[0.5 0.80 0.5],'edgecolor','none');

   % Biscayne National Park:
   m_line(bnp(:,1),bnp(:,2),'LineWidth',0.5,'LineStyle','-','Color',[1 1 1]*0.85);

   % Mangrove:
   m_bubble([LK.lon;ML.lon],[LK.lat;ML.lat],[meanLK{i};meanML{i}],mx*scaleM(i),1,1);

   % Reef (neg = black symbols):
   m_bubble(reef.lon,reef.lat,meanReef{i}*-1,mx*scaleR(i),1,1);

   % Annotate the map:
   m_text(-80.139,25.2575,{stage_txt{i}},'HorizontalAlignment','center',...
      'VerticalAlignment','bottom', 'FontSize',fnt,'Color','k');


   % -----Complete Base Map:-----
   m_grid('box','fancy','fontsize',8,'linestyle','none','xticklabels',[],...
      'xtick',(-80.375:[7.5/60]:-80.0625),'yticklabels',[],...
      'ytick',(25.25:[7.5/60]:25.6875));

end
