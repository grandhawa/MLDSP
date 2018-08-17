% Example of DISPROF based Cluster Analysis
% by David L. Jones, Jan-2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Fisher's Iris Data:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clz;

% Load the data::
load iris.mat;

% Perform cluster analysis via DISPROF + UPGMA:
dpa = f_disprof_clust(iris,'euc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 4566.2283, p = 0.0010
%   3 groups: Pi = 969.9449, p = 0.0010
%   4 groups: Pi = 150.2510, p = 0.0010
%   5 groups: Pi = 64.0851, p = 0.0010
%   6 groups: Pi = 59.4293, p = 0.0010
%   7 groups: Pi = 68.4920, p = 0.0010
%   8 groups: Pi = 2.2974, p = 0.2310
%   8 groups: Pi = 54.6769, p = 0.0010
%   9 groups: Pi = 8.8069, p = 0.0060
%   10 groups: Pi = 6.8857, p = 0.0270
%   11 groups: Pi = 11.4053, p = 0.1450
%   11 groups: Pi = 5.1859, p = 0.0320
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 33.5384, p = 0.0010
%   13 groups: Pi = 4.1157, p = 0.0290
%   14 groups: Pi = 4.2591, p = 0.0430
%   15 groups: Pi = 1.7625, p = 0.9670
%   15 groups: Pi = 0.7059, p = 0.3480
%   15 groups: Pi = 0.1983, p = 0.6170
%   15 groups: Pi = 2.7840, p = 0.0930
%   15 groups: Pi = 3.9998, p = 0.1570
%   15 groups: Pi = 1.7078, p = 0.0440
%   16 groups: Pi = 0.4099, p = 0.0270
%   17 groups: Pi = 0.1457, p = 0.6250
%   17 groups: Pi = 0.4924, p = 0.6420
%   17 groups: Pi = 0.1089, p = 0.2050
% 
% - # groups identified = 16 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = euc 
% -------------------------------------------------------
% 
% -> Sixteen groups in Fisher's Iris data? That's a lot of groups!

% Plot the cluster analysis tree as a dendrogram; branches depicted in gray
% indicate groups of observations that are homogeneous in terms of their
% composition and should not be split, but rather collapsed into a single
% cluster:
f_disprof_clustPlot(dpa,[],0,0.25);
% 
% -> note top=0 is used to create a dendrogram that spans right-to-left rather
% than top-to-bottom; also scale=0.25 is used to reduce the font size to
% eliminate overlap of the tic labels.

% Save as PDF:
f_pdf('iris_dendrogram');



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
% This example follows the analysis of Clarke et al. (2008) depicted in their
% Figure 4, but employs an equivalent approach based on disimilarity profile
% analysis (DISPROF) rather than similarity profile analysis (SIMPROF): 
% 
% Clarke, K. R., P. J. Somerfield, and R. N. Gorley. 2008.  Testing null
% hypotheses in exploratory community analyses: similarity profiles and
% biota-environmental linkage. J. Exp. Mar. Biol. Ecol. 366:56-69.

% Clear the workspace:
clz;

% Load the data::
load exe_nematodes.mat;

% 4th-root transform the data:
Y = f_normal(bio.dat,'4');

% Perform cluster analysis via DISPROF + UPGMA:
dpa = f_disprof_clust(Y,'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 17.1517, p = 0.0010
%   3 groups: Pi = 4.2715, p = 0.0010
%   4 groups: Pi = 4.8542, p = 0.0010
%   5 groups: Pi = 2.2085, p = 0.0010
%   6 groups: Pi = 1.8082, p = 0.0010
%   7 groups: Pi = 3.4130, p = 0.0010
%   8 groups: Pi = 0.6845, p = 0.0010
%   9 groups: Pi = 0.2014, p = 0.0300
%   10 groups: Pi = 0.0000, p = 1.0000
%   10 groups: Pi = 0.2089, p = 0.0010
%   11 groups: Pi = 0.0000, p = 1.0000
%   11 groups: Pi = 0.0954, p = 0.0480
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.1173, p = 0.0890
%   12 groups: Pi = 0.0000, p = 1.0000
% 
% - # groups identified = 11 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------

% Show resulting grouping vector:
dpa.grp
% ans =
% 
%      7
%      7
%      7
%      7
%      3
%      2
%     11
%     11
%      4
%      3
%      2
%     10
%     10
%      8
%      5
%      9
%      6
%      6
%      1

% Show incremental changes to grouping vector as branches in the cluster
% analysis tree are traversed by the DISPROF algorithm.
dpa.inc
% ans =
% 
%      1     2     2     4     4     4     7     7     7     7     7     7     7     7     7     7     7
%      1     2     2     4     4     4     7     7     7     7     7     7     7     7     7     7     7
%      1     2     2     4     4     4     7     7     7     7     7     7     7     7     7     7     7
%      1     2     2     4     4     4     7     7     7     7     7     7     7     7     7     7     7
%      1     1     3     3     3     3     3     3     3     3     3     3     3     3     3     3     3
%      1     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2
%      1     2     2     4     4     4     4     4     4     4     4     4    11    11    11    11    11
%      1     2     2     4     4     4     4     4     4     4     4     4    11    11    11    11    11
%      1     2     2     4     4     4     4     4     4     4     4     4     4     4     4     4     4
%      1     1     3     3     3     3     3     3     3     3     3     3     3     3     3     3     3
%      1     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2
%      1     1     1     1     1     1     1     8     8     8    10    10    10    10    10    10    10
%      1     1     1     1     1     1     1     8     8     8    10    10    10    10    10    10    10
%      1     1     1     1     1     1     1     8     8     8     8     8     8     8     8     8     8
%      1     1     1     1     5     5     5     5     5     5     5     5     5     5     5     5     5
%      1     1     1     1     1     6     6     6     9     9     9     9     9     9     9     9     9
%      1     1     1     1     1     6     6     6     6     6     6     6     6     6     6     6     6
%      1     1     1     1     1     6     6     6     6     6     6     6     6     6     6     6     6
%      1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1

% Show the rows of Y that were used in the DISPROF analysis to assess whether
% there were, say, 6 clusters in the data. Look at the above table and note that
% ' Evaluating the presence of...6 groups:' corresponds to the 5th DISPROF test
% performed, so show the rows associated with that assessment:
dpa.Pi(5)
dpa.p_idx{5}
% ans =
% 
%        1.8082
% 
% ans =
% 
%     12
%     13
%     14
%     16
%     17
%     18
%     19

% Plot the cluster analysis tree as a dendrogram; branches depicted in gray
% indicate groups of observations that are homogeneous in terms of their
% composition and should not be split, but rather collapsed into a single
% cluster:
f_disprof_clustPlot(dpa,[],1,1);

% Save as PDF:
f_pdf('exe_zooplankton_dendrogram');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Bristol Channel Plankton:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'bristol_zooplankton.mat' contains data from Collins & Williams
% (1982) concerning abundances (# per cubic meter) of 24 species of
% holozooplankton species collected from double oblique plankton net hauls at 57
% sites in the Bristol Channel, UK.
% 
% bio.dat   = abundance 24 spp of zooplankton from 57 sites
% bio.txt   = cell array of corresponding column (species) labels
% bio.sites = cell array of site labels
% 
% Collins, N. R. and R. Williams. 1982. Zooplankton communities in the Bristol
% Channel and Severn Estuary. Mar. Ecol. Prog. Ser. 9:1-11.
% 
% This example follows the analysis of Clarke et al. (2008) depicted in their
% Figure 4, but employs an equivalent approach based on disimilarity profile
% analysis (DISPROF) rather than similarity profile analysis (SIMPROF): 

% Clear the workspace:
clz;

% Load the data:
load bristol_zooplankton.mat;

% 4th-root transform the data:
Y = f_normal(bio.dat,'4');

% Perform cluster analysis via DISPROF:
dpa = f_disprof_clust(Y,'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 102.6760, p = 0.0010
%   3 groups: Pi = 35.9829, p = 0.0010
%   4 groups: Pi = 11.0068, p = 0.0010
%   5 groups: Pi = 1.0382, p = 0.0640
%   5 groups: Pi = 1.5334, p = 0.2630
%   5 groups: Pi = 0.7879, p = 0.3320
%   5 groups: Pi = 0.9488, p = 0.2450
% 
% - # groups identified = 4 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------

% Plot dendrogram:
f_disprof_clustPlot(dpa,bio.sites,0,0.65);
% -> note scale=0.65 is used to reduce the font size to eliminate overlap of the
% tic labels.

% Save as PDF:
f_pdf('bristol_zooplankton_dendrogram')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Ekofisk Oilfield Macrofauna:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'ekofisk_macrofauna.mat' consists of North Sea Ekofisk oil field
% data from Gray et al. (1988) and has the following variables: 
% 
% bio.dat   = abundance of 174 spp. of soft-bottom benthic macrofauna from 39 sites
% bio.txt   = cell array of corresponding column labels
% bio.sites = cell array of site labels
% 
% Gray, J. S., M. Aschan, M. R. Carr,  K. R. Clarke, R. H. Green, T. H. Pearson,
% R. Rosenberg, & R. M. Warwick. 1988. Analysis of community attributes of the
% benthic macrofauna of Frierfjord/Langesundfjord and in a mesocosm experiment.
% Mar. Ecol. Prog. Ser. 46: 151-165.
% 
% This example follows the analysis of Clarke et al. (2008) depicted in their
% Figure 4, but employs an equivalent approach based on disimilarity profile
% analysis (DISPROF) rather than similarity profile analysis (SIMPROF):

% Clear the workspace:
clz;

% Load the data:
load ekofisk_macrofauna.mat;

% Square-root transform the data:
Y = f_normal(bio.dat,'2');

% Perform cluster analysis via DISPROF:
dpa = f_disprof_clust(Y,'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of...
%   2 groups: Pi = 36.7424, p = 0.0010
%   3 groups: Pi = 0.0000, p = 1.0000
%   3 groups: Pi = 22.1904, p = 0.0010
%   4 groups: Pi = 8.1721, p = 0.0010
%   5 groups: Pi = 0.3512, p = 0.0840
%   5 groups: Pi = 2.2635, p = 0.0010
%   6 groups: Pi = 1.4711, p = 0.0010
%   7 groups: Pi = 0.1809, p = 0.0250
%   8 groups: Pi = 0.9608, p = 0.0090
%   9 groups: Pi = 0.8016, p = 0.0120
%   10 groups: Pi = 0.6763, p = 0.0260
%   11 groups: Pi = 0.3424, p = 0.0330
%   12 groups: Pi = 0.0000, p = 1.0000
%   12 groups: Pi = 0.0972, p = 0.1290
%   12 groups: Pi = 0.1895, p = 0.2850
%   12 groups: Pi = 0.0566, p = 0.2180
%   12 groups: Pi = 0.0000, p = 1.0000
% 
% - # groups identified = 11 
% - # permutation iterations = 1000 
% - alpha = 0.0500 
% - p-values adjusted for multiple testing: NO 
% - dissimilarity metric = bc 
% -------------------------------------------------------

% Plot dendrogram:
f_disprof_clustPlot(dpa,bio.sites,0,0.75);
% 
% -> note scale=0.75 is used to reduce the font size to eliminate overlap of the
% tic labels.

% Save as PDF:
f_pdf('ekofisk_macrofauna_dendrogram')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             NW Groundfish:                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clz;

% Load the data:
load nw_groundfish.mat;

% Square-root transform the data:
Y = f_normal(bio.dat,'2');

% Perform cluster analysis via DISPROF:
dpa = f_disprof_clust(Y,'bc');
% 
% =======================================================
%             Cluster Analysis via DISPROF:  
% -------------------------------------------------------
% Evaluating the presence of 2 groups...Pi = 1966.9225, p = 0.0010
% Evaluating the presence of 3 groups...Pi = 1905.0664, p = 0.0010
% Evaluating the presence of 4 groups...Pi = 1717.0121, p = 0.0010
% Evaluating the presence of 5 groups...Pi = 23.5410, p = 0.0010
% Evaluating the presence of 6 groups...Pi = 2.0408, p = 0.0010
% Evaluating the presence of 7 groups...Pi = 3.0253, p = 0.0010
% Evaluating the presence of 8 groups...Pi = 991.5312, p = 0.0010
% Evaluating the presence of 9 groups...Pi = 693.3954, p = 0.0010
% Evaluating the presence of 10 groups...Pi = 0.6791, p = 0.0010
% Evaluating the presence of 11 groups...Pi = 663.3052, p = 0.0010
% Evaluating the presence of 12 groups...Pi = 7.8226, p = 0.0010
% Evaluating the presence of 13 groups...Pi = 0.4383, p = 0.0390
% Evaluating the presence of 14 groups...Pi = 7.0092, p = 0.0010
% Evaluating the presence of 15 groups...Pi = 14.4216, p = 0.0010
% Evaluating the presence of 16 groups...Pi = 462.0909, p = 0.0010
% Evaluating the presence of 17 groups...Pi = 5.3854, p = 0.0010
% Evaluating the presence of 18 groups...Pi = 4.6517, p = 0.0010
% Evaluating the presence of 19 groups...Pi = 1.9097, p = 0.0010
% Evaluating the presence of 20 groups...Pi = 0.4709, p = 0.0020
% Evaluating the presence of 21 groups...Pi = 4.0509, p = 0.0010
% Evaluating the presence of 22 groups...Pi = 0.2872, p = 0.1380
% Evaluating the presence of 22 groups...Pi = 4.7280, p = 0.0010
% Evaluating the presence of 23 groups...Pi = 395.1750, p = 0.0010
% Evaluating the presence of 24 groups...Pi = 349.8759, p = 0.0010
% Evaluating the presence of 25 groups...Pi = 1.8513, p = 0.0010
% Evaluating the presence of 26 groups...Pi = 1.0203, p = 0.0040
% Evaluating the presence of 27 groups...Pi = 4.6922, p = 0.0010
% Evaluating the presence of 28 groups...Pi = 0.6641, p = 0.0030
% Evaluating the presence of 29 groups...Pi = 0.4674, p = 0.0050
% Evaluating the presence of 30 groups...Pi = 1.5616, p = 0.0020
% Evaluating the presence of 31 groups...Pi = 0.2664, p = 0.0050
% Evaluating the presence of 32 groups...Pi = 11.5609, p = 0.0010
% Evaluating the presence of 33 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 33 groups...Pi = 1.5502, p = 0.0010
% Evaluating the presence of 34 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 34 groups...Pi = 0.0555, p = 0.5510
% Evaluating the presence of 34 groups...Pi = 10.3778, p = 0.0010
% Evaluating the presence of 35 groups...Pi = 0.1284, p = 0.4280
% Evaluating the presence of 35 groups...Pi = 3.1612, p = 0.0020
% Evaluating the presence of 36 groups...Pi = 175.4030, p = 0.0010
% Evaluating the presence of 37 groups...Pi = 0.1177, p = 0.6850
% Evaluating the presence of 37 groups...Pi = 0.4923, p = 0.2000
% Evaluating the presence of 37 groups...Pi = 0.1321, p = 0.2350
% Evaluating the presence of 37 groups...Pi = 2.5240, p = 0.0010
% Evaluating the presence of 38 groups...Pi = 0.1964, p = 0.0020
% Evaluating the presence of 39 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 39 groups...Pi = 1.6677, p = 0.0480
% Evaluating the presence of 40 groups...Pi = 0.5831, p = 0.0510
% Evaluating the presence of 41 groups...Pi = 2.9278, p = 0.0010
% Evaluating the presence of 42 groups...Pi = 1.2125, p = 0.0680
% Evaluating the presence of 42 groups...Pi = 0.3148, p = 0.0130
% Evaluating the presence of 43 groups...Pi = 68.8816, p = 0.0010
% Evaluating the presence of 44 groups...Pi = 0.1245, p = 0.9250
% Evaluating the presence of 44 groups...Pi = 0.5570, p = 0.0430
% Evaluating the presence of 45 groups...Pi = 0.0665, p = 0.2760
% Evaluating the presence of 45 groups...Pi = 0.2251, p = 0.0390
% Evaluating the presence of 46 groups...Pi = 0.1322, p = 0.0660
% Evaluating the presence of 46 groups...Pi = 18.2188, p = 0.0010
% Evaluating the presence of 47 groups...Pi = 0.0905, p = 0.1020
% Evaluating the presence of 47 groups...Pi = 14.4211, p = 0.0010
% Evaluating the presence of 48 groups...Pi = 1.2789, p = 0.0050
% Evaluating the presence of 49 groups...Pi = 0.1145, p = 0.0540
% Evaluating the presence of 50 groups...Pi = 0.1720, p = 0.1280
% Evaluating the presence of 50 groups...Pi = 0.8099, p = 0.0070
% Evaluating the presence of 51 groups...Pi = 10.8553, p = 0.0010
% Evaluating the presence of 52 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 52 groups...Pi = 9.7546, p = 0.0010
% Evaluating the presence of 53 groups...Pi = 0.0237, p = 0.8680
% Evaluating the presence of 53 groups...Pi = 0.2613, p = 0.2470
% Evaluating the presence of 53 groups...Pi = 0.3498, p = 0.0650
% Evaluating the presence of 53 groups...Pi = 0.4403, p = 0.0110
% Evaluating the presence of 54 groups...Pi = 0.2319, p = 0.6240
% Evaluating the presence of 54 groups...Pi = 0.1614, p = 0.2240
% Evaluating the presence of 54 groups...Pi = 6.0919, p = 0.0010
% Evaluating the presence of 55 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 55 groups...Pi = 0.1300, p = 0.1850
% Evaluating the presence of 55 groups...Pi = 1.4611, p = 0.0010
% Evaluating the presence of 56 groups...Pi = 0.3361, p = 0.0340
% Evaluating the presence of 57 groups...Pi = 0.0306, p = 0.6410
% Evaluating the presence of 57 groups...Pi = 0.5663, p = 0.0040
% Evaluating the presence of 58 groups...Pi = 0.0217, p = 0.8880
% Evaluating the presence of 58 groups...Pi = 0.0893, p = 0.5810
% Evaluating the presence of 58 groups...Pi = 0.2475, p = 0.0450
% Evaluating the presence of 59 groups...Pi = 0.0430, p = 0.4780
% Evaluating the presence of 59 groups...Pi = 4.1666, p = 0.0010
% Evaluating the presence of 60 groups...Pi = 0.2347, p = 0.0480
% Evaluating the presence of 61 groups...Pi = 3.7427, p = 0.0010
% Evaluating the presence of 62 groups...Pi = 0.1374, p = 0.4410
% Evaluating the presence of 62 groups...Pi = 7.3392, p = 0.0010
% Evaluating the presence of 63 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 63 groups...Pi = 0.2683, p = 0.0060
% Evaluating the presence of 64 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 64 groups...Pi = 0.0230, p = 0.7820
% Evaluating the presence of 64 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 64 groups...Pi = 0.0536, p = 0.5770
% Evaluating the presence of 64 groups...Pi = 0.4502, p = 0.1470
% Evaluating the presence of 64 groups...Pi = 0.0420, p = 0.5510
% Evaluating the presence of 64 groups...Pi = 0.0531, p = 0.2820
% Evaluating the presence of 64 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 64 groups...Pi = 0.0409, p = 0.4270
% Evaluating the presence of 64 groups...Pi = 0.8821, p = 0.0070
% Evaluating the presence of 65 groups...Pi = 1.3497, p = 0.0020
% Evaluating the presence of 66 groups...Pi = 0.3774, p = 0.0020
% Evaluating the presence of 67 groups...Pi = 0.1565, p = 0.0480
% Evaluating the presence of 68 groups...Pi = 6.4953, p = 0.0010
% Evaluating the presence of 69 groups...Pi = 5.1872, p = 0.0010
% Evaluating the presence of 70 groups...Pi = 0.0540, p = 0.7530
% Evaluating the presence of 70 groups...Pi = 0.8174, p = 0.0080
% Evaluating the presence of 71 groups...Pi = 0.1536, p = 0.1910
% Evaluating the presence of 71 groups...Pi = 0.9337, p = 0.0020
% Evaluating the presence of 72 groups...Pi = 0.2767, p = 0.0790
% Evaluating the presence of 72 groups...Pi = 0.4518, p = 0.0460
% Evaluating the presence of 73 groups...Pi = 0.0608, p = 0.0610
% Evaluating the presence of 73 groups...Pi = 0.0451, p = 0.2290
% Evaluating the presence of 73 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 73 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 73 groups...Pi = 0.5906, p = 0.0050
% Evaluating the presence of 74 groups...Pi = 0.0280, p = 0.6790
% Evaluating the presence of 74 groups...Pi = 0.2779, p = 0.0470
% Evaluating the presence of 75 groups...Pi = 0.4647, p = 0.0030
% Evaluating the presence of 76 groups...Pi = 0.0774, p = 0.7020
% Evaluating the presence of 76 groups...Pi = 0.0517, p = 0.0690
% Evaluating the presence of 76 groups...Pi = 0.0648, p = 0.0480
% Evaluating the presence of 77 groups...Pi = 0.0401, p = 0.2410
% Evaluating the presence of 77 groups...Pi = 0.1462, p = 0.1700
% Evaluating the presence of 77 groups...Pi = 0.0404, p = 0.8280
% Evaluating the presence of 77 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 77 groups...Pi = 0.2809, p = 0.6390
% Evaluating the presence of 77 groups...Pi = 0.2488, p = 0.1220
% Evaluating the presence of 77 groups...Pi = 0.0519, p = 0.3580
% Evaluating the presence of 77 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 77 groups...Pi = 0.0544, p = 0.2440
% Evaluating the presence of 77 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 77 groups...Pi = 0.0000, p = 1.0000
% Evaluating the presence of 77 groups...Pi = 0.0000, p = 1.0000
% 
% Note: 76 natural groups have been identified in the data
% using 1000 permutation iterations and alpha = 0.0500.
% (p-values were NOT corrected for multiple testing)
% -------------------------------------------------------

% Plot dendrogram:
f_disprof_clustPlot(dpa,bio.sites,0,1/8);
% 
% -> note top=0 is used to create a dendrogram that spans right-to-left rather
% than top-to-bottom; also scale=1/8 is used to dramatically reduce the font
% size to eliminate overlap of the tic labels.

% Save as PDF:
f_pdf('nw_groundfish_dendrogram')
