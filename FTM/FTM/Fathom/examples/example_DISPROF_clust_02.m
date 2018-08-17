% More examples of DISPROF based Cluster Analysis
% 
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          EXE Estuary Nematodes:                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'exe_nematodes.mat' consists of abundance data for 140 species of
% nematodes from 19 sites within the Exe estuary in the UK from Warwick (1971)
% and has the following variables:
% 
% bio.dat   = average counts of 140 spp of nematodes from 19 sites
% bio.txt   = cell array of corresponding column (species) labels
% bio.sites = cell array of site labels
% 
% Warwick, R. M. 1971. Nematode associations in the Exe estuary. J. Mar. Biol.
% Assoc. U.K. 51: 439-454.
% 

% Clear the workspace:
clz;

% Load the data::
load exe_nematodes.mat;

% 4th-root transform the data:
Y = f_normal(bio.dat,'4');

% Perform cluster analysis via DISPROF + UPGMA + Bray-Curtis metric:
[dpa_1,bin_1] = f_disprof_clust(Y,'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 17.1271, p = 0.0010
%   3 groups: Pi = 4.2954, p = 0.0010
%   4 groups: Pi = 4.8487, p = 0.0010
%   5 groups: Pi = 2.2003, p = 0.0010
%   6 groups: Pi = 1.7917, p = 0.0010
%   7 groups: Pi = 3.4064, p = 0.0010
%   8 groups: Pi = 0.6859, p = 0.0010
%   9 groups: Pi = 0.1978, p = 0.0350
%   10 groups: Pi = 0.0000, p = 1.0000
%   10 groups: Pi = 0.2079, p = 0.0010
%   11 groups: Pi = 0.0000, p = 1.0000
%   11 groups: Pi = 0.0990, p = 0.0510
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.1189, p = 0.0990
%   12 groups: Pi = 0.0000, p = 1.0000
% 
% - # groups identified = 11 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------


% Perform cluster analysis via DISPROF + WPGMA + Bray-Curtis metric:
[dpa_2,bin_2] = f_disprof_clust(Y,'bc',6);
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 17.1394, p = 0.0010
%   3 groups: Pi = 7.4900, p = 0.0010
%   4 groups: Pi = 4.8528, p = 0.0010
%   5 groups: Pi = 2.1989, p = 0.0010
%   6 groups: Pi = 0.4607, p = 0.0130
%   7 groups: Pi = 3.4139, p = 0.0010
%   8 groups: Pi = 0.6846, p = 0.0010
%   9 groups: Pi = 0.1982, p = 0.0420
%   10 groups: Pi = 0.0000, p = 1.0000
%   10 groups: Pi = 0.2076, p = 0.0010
%   11 groups: Pi = 0.0000, p = 1.0000
%   11 groups: Pi = 0.0950, p = 0.0480
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.1221, p = 0.0680
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.0000, p = 1.0000
% 
% - # groups identified = 11 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------


% Compare the dendrograms by the 2 different clustering methods:
f_disprof_clustPlot(dpa_1,[],0);
f_disprof_clustPlot(dpa_2,[],0);
% 
% -> they look pretty similar, lets check the grouping vectors:
[(1:numel(dpa_1.grp))' dpa_1.grp dpa_2.grp]
% ans =
% 
%      1     7     7
%      2     7     7
%      3     7     7
%      4     7     7
%      5     3     3
%      6     2     2
%      7    11    11
%      8    11    11
%      9     4     4
%     10     3     3
%     11     2     2
%     12    10    10
%     13    10    10
%     14     8     8
%     15     5     6
%     16     9     9
%     17     6     5
%     18     6     5
%     19     1     1
% 
% -> Some of the cluster designations are different: object 15 was assigned to
% cluster 5 using UPGMA but to cluster 6 using WPGMA. Also, object 17 &
% 18 were assigned to cluster 6 using UPGMA, but to cluster 5 using WPGMA
% method. However the overall structure remains the same, only the cluster 
% 'name' has changed.

% Measure the congruence between the two clustering methods by comparing the
% binary connectivity matrix prodced by each:
U_W = f_disprof_clust_bin(bin_1,bin_2)
% U_W = 
% 
%       n: 0
%     idx: []
%     cng: 1
% 
% -> congruence ranges from 0-1, so in this case there is 100% agreement between
%    the two methods; any differences in the way the clusters were 'named'
%    between the two methods arose due to the ORDER clusters were created and
%    not in their underlying membership.


% Perform cluster analysis via DISPROF + UPGMA + Euclidean Distance:
[dpa_3,bin_3] = f_disprof_clust(Y,'euc',1);
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 240.0329, p = 0.0010
%   3 groups: Pi = 202.7713, p = 0.0010
%   4 groups: Pi = 122.5597, p = 0.0010
%   5 groups: Pi = 71.6664, p = 0.0010
%   6 groups: Pi = 1.3642, p = 0.0800
%   6 groups: Pi = 3.4294, p = 0.0010
%   7 groups: Pi = 25.6985, p = 0.0010
%   8 groups: Pi = 13.4229, p = 0.0010
%   9 groups: Pi = 5.0276, p = 0.0030
%   10 groups: Pi = 0.0000, p = 1.0000
%   10 groups: Pi = 1.2444, p = 0.0960
%   10 groups: Pi = 1.1959, p = 0.1530
%   10 groups: Pi = 0.0000, p = 1.0000
% 
% - # groups identified = 9 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = euc 
% -------------------------------------------------------

% View dendrogram:
f_disprof_clustPlot(dpa_3,[],0);


% Compare grouping vectors:
[(1:numel(dpa_1.grp))' dpa_1.grp dpa_3.grp]
% 
% ans =
% 
%      1     7     5
%      2     7     5
%      3     7     5
%      4     7     5
%      5     3     9
%      6     2     7
%      7    11     3
%      8    11     3
%      9     4     3
%     10     3     9
%     11     2     2
%     12    10     6
%     13    10     6
%     14     8     4
%     15     5     8
%     16     9     8
%     17     6     8
%     18     6     8
%     19     1     1

% Measure the congruence between the two clustering methods by comparing the
% binary connectivity matrix prodced by each:
bc_euc = f_disprof_clust_bin(bin_1,bin_3)
% bc_euc = 
%       n: 5
%     idx: [5x1 double]
%     cng: 0.73684
% 
% -> 5 of the 19 objects in Y are assigned to different clusters by these 2
%    methods, so the level of agreement is 74%
