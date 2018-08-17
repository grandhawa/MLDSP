% Example of Disimilarity Profile Analysis (DISPROF, SIMPROF)
% 
% by David L. Jones, Jan-2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Exe Estuary Nematodes:                           %
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

% This example follows the analysis of Clarke et al. (2008) depicted in their
% Figure 5 but employs an equivalent approach based on disimilarity profile
% analysis (DISPROF) rather than similarity profile analysis (SIMPROF). DISPROF
% is used to test the null hypothesis that there is no multivariate structure in
% the abundance and distribution data represented by the Exe data set.

% Clear the workspace:
clz;

% Load data:
load exe_nematodes.mat;

% Fourth-root transform the abundance data:
bio.dat4 = f_normal(bio.dat,'4');

% Run DISPROF using the Bray-Curtis disimilarity metric:
dpa = f_disprof(bio.dat4,'bc',1000,1,1);
% 
% Calculating mean disimilarity profile from 999 permutations...
% Calculating significance levels from 999 permutations...
% 
% ==================================================
%       DISPROF - Disimilarity Profile Analysis:
% --------------------------------------------------
% Pi stat = 17.1427  p =  0.00100 
% No. of permutations = 1000 
% ==================================================
% 
% -> the Pi statistic is highly significant, so we reject the null hypothesis of
%    no significant multivariate structure in the data (at alpha = 0.05). Note
%    the observed disimilarity profile depicted in the corresponding plot
%    departs substantially from the mean profile generated through random
%    permutation. The observed profile has more small and large disimilarity
%    values than that expected by chance, suggesting there is true strcture in
%    these data.

% Repeat the assessment, but only for sites 1-4:
dpa_14 = f_disprof(bio.dat4(1:4,:),'bc',1000,1,1);
% 
% Calculating mean disimilarity profile from 999 permutations...
% Calculating significance levels from 999 permutations...
% 
% ==================================================
%       DISPROF - Disimilarity Profile Analysis:
% --------------------------------------------------
% Pi stat = 0.1182  p =  0.08800 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% -> By limiting the analysis to just sites 1-4, the Pi statistic is not
% significant, so we accept the null hypothesis of 'Homogeneity of Composition'
% among these sites (at alpha = 0.05) and consider this group of sites to be
% members of single, homogeneous cluster.

