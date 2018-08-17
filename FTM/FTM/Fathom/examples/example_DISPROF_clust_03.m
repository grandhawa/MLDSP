% More examples of DISPROF based Cluster Analysis
% 
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tikus Coral Assemblages:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'tikus.mat' contains data on coral assemblages from Tikus
% Island, Indonesia taken along 10 replicate 30 m line transects during
% each of 6 years (Warwick et al., 1990). The 1982 El Niño triggered a
% disturbance event which resulted in widespread coral bleaching in this
% region. As such, the 1981, 1983, and 1985 samples were collected before,
% during, and after the El Niño disturbance event, repectively.
% 
% Warwick, R. M., K. R. Clarke, and Suharsono. 1990. A statistical analysis
% of coral community responses to the 1982-1983 El Niño in the Thousand
% Islands, Indonesia. Coral Reefs 8:171-179.
% 
% -----Variables:----
% coral     = percent cover of 75 species in each transect
% coral_txt = cell array of corresponding column labels
% yr        = year (1 = 1981, 3 = 1983, 4 = 1984, 5 = 1985, 7 = 1987, 8 = 1988)
% yr_txt    = cell array of corresponding row labels

% 1) Examine the effect of this disturbance event on the coral assemblages
% by determining whether there were differences in the variation in the
% abundance and composition of corals before (1981), during (1983), and
% after (1985) the 1982 El Niño event. Interpret your results in terms of a
% null hypothesis.

% Load data:
load tikus.mat

% Apply a mild data transformation to down-weight the influence of the most
% abundant species:
coral_2 = f_normal(coral,'2');

% Perform cluster analysis via DISPROF + UPGMA + Bray-Curtis metric:
[dpa_1,bin_1] = f_disprof_clust(coral_2,'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 73.7066, p = 0.0010
%   3 groups: Pi = 65.2450, p = 0.0010
%   4 groups: Pi = 0.4829, p = 0.4660
%   4 groups: Pi = 49.0355, p = 0.0010
%   5 groups: Pi = 27.6831, p = 0.0010
%   6 groups: Pi = 6.3471, p = 0.0010
%   7 groups: Pi = 5.6338, p = 0.0290
%   8 groups: Pi = 0.4004, p = 0.2570
%   8 groups: Pi = 1.8959, p = 0.0010
%   9 groups: Pi = 4.9283, p = 0.0360
%   10 groups: Pi = 3.1145, p = 0.0030
%   11 groups: Pi = 4.8763, p = 0.0030
%   12 groups: Pi = 1.2005, p = 0.6930
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 1.3324, p = 0.0020
%   13 groups: Pi = 0.9689, p = 0.0030
%   14 groups: Pi = 0.1373, p = 0.1220
%   14 groups: Pi = 0.1939, p = 0.3040
%   14 groups: Pi = 0.1863, p = 0.1960
%   14 groups: Pi = 0.2124, p = 0.2030
%   14 groups: Pi = 0.0000, p = 1.0000
%   14 groups: Pi = 0.2545, p = 0.2930
%   14 groups: Pi = 0.1549, p = 0.4890
% 
% - # groups identified = 13 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------

% View the dendrogram:
f_disprof_clustPlot(dpa_1,yr_txt,0);


% Re-do the analysis, but exclude data collected during the El Niño event:
idx           = yr~=3; % get index to samples not collected during 1983s
[dpa_2,bin_2] = f_disprof_clust(coral_2(idx,:),'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 50.7263, p = 0.0010
%   3 groups: Pi = 47.1339, p = 0.0010
%   4 groups: Pi = 26.4876, p = 0.0010
%   5 groups: Pi = 4.9745, p = 0.0010
%   6 groups: Pi = 1.9063, p = 0.0010
%   7 groups: Pi = 3.1471, p = 0.0010
%   8 groups: Pi = 5.5408, p = 0.0060
%   9 groups: Pi = 0.1538, p = 0.6040
%   9 groups: Pi = 0.4721, p = 0.2870
%   9 groups: Pi = 1.3184, p = 0.6720
%   9 groups: Pi = 1.3447, p = 0.0010
%   10 groups: Pi = 0.9694, p = 0.0030
%   11 groups: Pi = 0.0000, p = 1.0000
%   11 groups: Pi = 0.1399, p = 0.1250
%   11 groups: Pi = 0.2194, p = 0.2080
%   11 groups: Pi = 0.0000, p = 1.0000
%   11 groups: Pi = 0.2513, p = 0.3090
%   11 groups: Pi = 0.1498, p = 0.5030
% 
% - # groups identified = 10 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------

% View the dendrogram:
f_disprof_clustPlot(dpa_2,yr_txt(idx),0);



% See how much of an effect removal of the 1983 samples had on clustering of the
% remaining samples. Do this by measuring the congruence between the results of
% the two cluster analyses by comparing the binary connectivity matrices (BIN)
% prodced by each.
% 
% First, generate an index that specifies the rows both analyses have in common:
idx_1 = find(ismember(yr,yr(idx)));

% Second, compare a subset of the first BIN matrix with the corresponding
% elements of the second BIN matrix:
result = f_disprof_clust_bin(bin_1(idx_1,idx_1),bin_2)
% result = 
% 
%       n: 20
%     idx: [20x1 double]
%     cng: 0.6

% -> 20 of the 50 objects were assigned to different clusters, so the level of
% agreement between the two cluster analyses is 60%.
