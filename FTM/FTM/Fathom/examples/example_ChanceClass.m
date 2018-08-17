% Examples of Classification Success Rates Expected by Chance
% 
% by David L. Jones, May-2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% What is the probability of correctly classifying a group of samples by
% chance alone? For equal sample sizes, the probability of correct
% classification is 1/g, where g = # groups
% 
y = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
f_chanceClass(y);
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1             6.25 % 
%    2             6.25 % 
%    3             6.25 % 
%    4             6.25 % 
% 
% Total Correct = 25.00 % 
% --------------------------------------------------

% For highly unbalanced data, where one group tends to dominate, the chance
% classification success rate increases towards 100%:
w = [1 2 3 ones(1,99)*4];
f_chanceClass(w);
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1             0.01 % 
%    2             0.01 % 
%    3             0.01 % 
%    4            94.20 % 
% 
% Total Correct = 94.23 % 
% --------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Galapagos fish:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Ruttenberg & Warner (2006) used LA-ICP-MS to measure the otolith
% microchemistry of 119 embryos from benthic egg clutches of damselfish
% from three regions in the western Galapagos islands. Their samples sizes
% were n = 50, 45, and 24 for each of the three regions, respectively. A
% discriminat analysis yielded a classifier that had a 66.4%
% re-classification success rate. Ruttenberg, B. I. and R. R. Warner. 2006.
% Variation in the chemical composition of natal otoliths of a reef fish
% from the Galapagos Islands. Mar. Ecol. Prog. Ser. 328: 225-236.
% 
% Was this better than that expected from chance alone?
% 
grp = [ones(50,1); ones(45,1)*2; ones(24,1)*3];
f_chanceClass(grp,0.664,1000,1);
% 
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1            17.65 % 
%    2            14.30 % 
%    3             4.07 % 
% 
% Total Correct = 36.02 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 35.84 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% -> reject the null hypothesis that the observed re-classification success
%    rate was no better than that expected by chance alone, at alpha =
%    0.05.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Simulated Otolith Microchemistry:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% White & Ruttenberg generated simulated otolith microchemistry data to
% examine the effect of relative sample sizes on classification accuracy.
% In one of their examples (Table 1), they generated data for 500
% individuals from 3 sites, with n = 450, 25, and 25 samples from each
% site, respectively. The re-classification success rate of the resulting
% discriminant analysis was 72.9%. White, J. W. and B. I. Ruttenberg. 2007.
% Discriminant function analysis in marine ecology: some oversights and
% their solutions. Mar. Ecol. Prog. Ser. 329: 301-305.
% 
% Was this better than that expected from chance alone?
grp = [ones(450,1); ones(25,1)*2; ones(25,1)*3];
f_chanceClass(grp,0.729,1000,1);
% 
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1            81.00 % 
%    2             0.25 % 
%    3             0.25 % 
% 
% Total Correct = 81.50 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 81.51 % 
% p =  1.00000 
% No. of permutations = 1000 
% --------------------------------------------------
% -> ACCEPT the null hypothesis that the observed re-classification success
%    rate was no better than that expected by chance alone, at alpha =
%    0.05.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Leatherback Turtles:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Solow (1990) analyzed data concerning the genetic differences among 2
% subpopulations of leatherback turtles based on blood proteins. Data were
% collected from 22 nesting females: 12 turtles in Costa Rica and 10 in the
% Virgin Islands. A discriminant analysis with backwards elimination
% yielded a classifier with a re-classification success rate of 68.2%.
% Solow, A. R. 1990. A randomization test for misclassification probability
% in discriminant analysis. Ecology 71: 2379-2382.
% 
% Was their classification rule significantly better than random
% assignment?
grp = [ones(12,1); ones(10,1)*2;];
f_chanceClass(grp,0.682,5000,1);
% 
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1            29.75 % 
%    2            20.66 % 
% 
% Total Correct = 50.41 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 50.49 % 
% p =  0.04780 
% No. of permutations = 5000 
% --------------------------------------------------
% -> The p-value was marginally significant at the 0.05 level, so we REJECT
%    the null hypothesis that the observed re-classification success rate was
%    no better than that expected by chance alone.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Rainbow Parrotfishes:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Machemer et al. (2012) examined mangrove visual survey data collected in
% Biscayne Bay from 1998-2009. Rainbow parrotfish were found present in 106
% of 1779 sites surveyed during 1998-2007. Using a number of environmental
% variables and logistic regression with AIC-based variable selection, they
% generated a model that correctly classified 94.5% of the 2008-2009 data, 
% where parrotfish were observed in 53 of 759 sites. Machemer, E. G. P.,
% J. F. Walter III, J. E. Serafy, and D. W. Kerstetter. 2012. Importance of
% mangrove shorelines for rainbow parrotfish Scaurs guacamaia: habitat
% suitability modeling in a subtropical bay. Aquat. Biol. 15: 87-98.
% 
% Is their model better than that expected from simple random assignment?
grp = [ones(53,1); zeros(706,1)];
f_chanceClass(grp,0.945,1000,1);
% 
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1             0.49 % 
%    0            86.52 % 
% 
% Total Correct = 87.01 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 87.01 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% -> The p-value was highly significant at the 0.05 level, so we reject
%    the null hypothesis that the observed classification success rate was
%    no better than that expected by chance alone. However, we note that
%    their model's overall success rate was only 7.5 percentage points
%    better than random allocation.
