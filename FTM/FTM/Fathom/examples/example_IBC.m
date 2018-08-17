% exampleIBC.m
% - example of Indicator-species Based Cluster (IBC) analysis
% 
% by David L. Jones, Jan-2013

% Mar-2013: updated

% -----Notes:-----
% This is data is from from van der Aart & Smeenk-Enserink (1975) and is Hunting
% spider abundances for 12 species (spiders) taken from 28 sites (site_labels)
% and associated environmental data (env). 
 
% Load data file: 
load spiders.mat

% Create dissimilarity matrix:
dis = f_dis(spiders,'bc');

% Obtain a Euclidean-embedding via PCoA:
pcoa = f_pcoa(dis,0,1,0);

% Plot the ordination:
[nul,vec] = f_pcoaPlot(pcoa,[],spiders,spiders_labels,0,'none',1);

% Examine the species biplot vector values:
vec(:,1:2)
% 
% ans =
% 
%       0.29882     -0.18835
%      -0.18586     -0.33227
%        0.3243     0.034857
%       0.18877     0.059771
%      -0.06466     -0.16284
%      -0.19505     -0.43571
%      -0.22907      0.23219
%      0.093162     -0.59803
%      -0.31976      -0.7388
%       -0.4032           -1
%       -0.6725     -0.57695
%      -0.28882     -0.37706


% -----One way to select the 3 most influential species:-----
% Get length of each vector:
eDis = f_dis([0 0;vec(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Sort species according to length of biplot vector:
[nul,idxVec] = sortrows(eDis,-1);
spiders_labels(idxVec(1:3))' % list species

% ans = 
% 
%     'Pard_pull'
%     'Troc_terr'
%     'Pard_nigr'
% ------------------------------------------------

% Create a Minimum Spanning Tree:
mst = f_mst(dis,pcoa.scores);

% Plot the MST with abundance labels:
txt = f_num2cell(spiders(:,idxVec(2))); % labels specifying abundance of Troc_terr
f_plotNeigh(mst,2,'k',txt);
title('\bfMinimum Spanning Tree');
xlabel('Axis I');
ylabel('Axis II');
axis([-0.4 0.7 -0.4 0.5]);
f_pdf('mst_abund'); % save pdf

% Plot the MST with site labels:
f_plotNeigh(mst,2,'k');
title('\bfMinimum Spanning Tree');
xlabel('Axis I');
ylabel('Axis II');
axis([-0.4 0.7 -0.4 0.5]);
f_pdf('mst_site'); % save pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     IBC based on 3 most influential taxa:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
ibc = f_ibc(mst.con,spiders(:,idxVec(1:3)),'fileOUT',spiders_labels(:,idxVec(1:3)),site_labels,[])
% 
% ibc = 
% 
%        Y: [28x3 double]
%      tol: [0.6229 0.65322 0.69386]
%      con: {[7x2 double]  [8x2 double]  [2x2 double]}
%     idxS: {[8x1 double]  [9x1 double]  [3x1 double]}
%        S: [28x3 double]
%        B: {[28x28 double]  [28x28 double]  [28x28 double]}
%     uGrp: [3x3 double]
%     yTxt: [1x78 char]
%     sTxt: [1x48 char]

% Show the sites indicative of 'Troc_terr':
site_labels(ibc.idxS{2})
% ans = 
% 
%     's03'
%     's13'
%     's07'
%     's14'
%     's01'
%     's05'
%     's04'
%     's06'
%     's02'
% 
% Show the corresponding Boolean matrix of sites defining the cluster:
ibc.S(:,2)
% 
% ans =
% 
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      0
%      0
%      0
%      0
%      0
%      1
%      1
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
% 
% Determine the consensus cluster of sites for all 3 indicator species:
find( (ibc.S(:,1) .* ibc.S(:,2) .* ibc.S(:,3)) == 1) % "lowest common denominator"
% ans =
% 
%      5
%      7
%     13
% 
find( (ibc.S(:,1) + ibc.S(:,2) + ibc.S(:,3)) >0) % "aggregation"
% ans =
% 
%      1
%      2
%      3
%      4
%      5
%      6
%      7
%     13
%     14
% 
% Show list of groups of indicator species exported to *.csv file:
type fileOUT_y.csv
% 
% Pard_pull, Troc_terr, -, 
% Pard_pull, Troc_terr, Pard_nigr, 
% -, Troc_terr, -, 
% -, -, -, 
% 
% Show corresponding list of sites where these groups occur exported to *.csv file:
type fileOUT_s.csv
% 
% s01, s02, s03, s04, s14, 
% s05, s07, s13, 
% s06, 
% s08, s09, s10, s11, s12, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%             IBC based on all taxa:           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
ibc_2 = f_ibc(mst.con,spiders(:,idxVec(1:end)),'fileOUT_2',spiders_labels(:,idxVec(1:end)),site_labels,[])
% 
% ibc_2 = 
% 
%        Y: [28x12 double]
%      tol: [0.6229 0.65322 0.69386 0.70603 0.67772 0.64521 0.74606 0.66667 0.65353 0.77647 0.76013 0.75633]
%      con: {1x12 cell}
%     idxS: {1x12 cell}
%        S: [28x12 double]
%        B: {1x12 cell}
%     uGrp: [10x12 double]
%     yTxt: [1x594 char]
%     sTxt: [1x95 char]
% 
% 
% 
% Show unique rows of 'S':
ibc_2.uGrp
% 
% ans =
% 
%      1     1     0     1     0     0     0     1     0     0     0     0
%      1     1     0     0     0     0     0     0     0     0     0     0
%      1     1     0     0     1     1     0     0     0     0     0     0
%      1     1     1     0     0     1     0     0     0     0     0     0
%      0     1     0     0     0     1     0     0     0     0     0     0
%      1     1     1     0     1     1     0     0     0     0     0     0
%      0     0     0     1     0     0     0     0     0     0     0     0
%      0     0     0     1     0     0     0     1     0     0     0     0
%      1     1     0     0     0     1     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     1     0     0     0
% 
% 
% Show the indicator species comprising each unique group of taxa:
ibc_2.yTxt
% 
% ans =
% 
% Pard_pull, Troc_terr, -, Pard_mont, -, -, -, Alop_acce, -, -, -, -, 
% Pard_pull, Troc_terr, -, -, -, -, -, -, -, -, -, -, 
% Pard_pull, Troc_terr, -, -, Aulo_albi, Zora_spin, -, -, -, -, -, -, 
% Pard_pull, Troc_terr, Pard_nigr, -, -, Zora_spin, -, -, -, -, -, -, 
% -, Troc_terr, -, -, -, Zora_spin, -, -, -, -, -, -, 
% Pard_pull, Troc_terr, Pard_nigr, -, Aulo_albi, Zora_spin, -, -, -, -, -, -, 
% -, -, -, Pard_mont, -, -, -, -, -, -, -, -, 
% -, -, -, Pard_mont, -, -, -, Alop_acce, -, -, -, -, 
% Pard_pull, Troc_terr, -, -, -, Zora_spin, -, -, -, -, -, -, 
% -, -, -, -, -, -, -, -, Alop_fabr, -, -, -, 
% 
% 
% Show the sites where each unique group occurs:
ibc_2.sTxt
% 
% ans =
% 
% s01, 
% s02, s03, 
% s04, 
% s05, s13, 
% s06, 
% s07, 
% s11, 
% s12, 
% s14, 
% s22, s23, s24, s25, s27, s28, 
