% Examples for Chi-Square Test of Independence
% (= Homogeneity of Proportions)
% 
% by David L. Jones, 2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 1) Test the null hypothesis that there is no difference in habitat usage
% between two species of eels. Data are from Young and Winn (2003) and
% consist of the number of times each species was found within one of three
% habitat types:
% 
%             G. moringa:     G. vicinus:
% Grass       127 (25.9%)     116 (33.7%) 
% Sand         99 (20.2%)      67 (19.5%)
% Border      264 (53.95)     161 (46.8%)

% Create RxC Data table:
Y = [127 116; 99 67; 264 161];

% Chi-Square Test:
f_chisqIndep(Y,1000)
% ans = 
% 
% ans = 
% 
%       stat: 6.2621
%         df: 2
%     p_perm: 0.056
%     p_para: 0.043671
% 
% -> We reject the null hypothesis of no difference in habitat usage (at
% alpha = 0.05)

% Re-do the analysis, this time using proportions:
Y_freq = (Y./repmat(sum(Y),size(Y,1),1)) % convert counts to proportions
% Y_freq =
% 
%       0.25918      0.33721
%       0.20204      0.19477
%       0.53878      0.46802
% 
f_chisqIndep(Y_freq,1000)
% 
% ans = 
% 
%       stat: 0.015313
%         df: 2
%     p_perm: 0.052
%     p_para: 0.99237
% 
% -> this doesn't work, input data must be COUNTS, not proportions!

% 2) Test the null hypothesis that there is no difference in substrate
% usage between herons and egrets. Data are from Custer and Galli (2002)
% and consist of the number of times each species was recorded to land in
% each of 4 substrate types:
% 
%              Heron   Egret
% Vegetation    15      8
% Shoreline     20      5
% Water         14      7
% Structures     6      1

% Create RxC Data table:
X = [15 20 14 6; 8 5 7 1];

%  Chi-Square Test:
f_chisqIndep(X,1000)
% 
% ans = 
% 
%       stat: 2.2812
%         df: 3
%     p_perm: 0.111
%     p_para: 0.51613


% Slippery dick              10             6
% Unidentified wrasses        3             7
% Moray eels                  1             1
% Squirrelfish                1             1
% Unidentified fish           6             3

% Oval urn crab              31            10
% Emerald crab                3             2
% Portunus crab spp.          1             0
% Arrow crab                  1             0
% Unidentified crabs         15             1
% Spiny lobster               0             1

% Octopus                     3             2
% Unidentified                4             1
% 
% 
% The nominal variables are the species of eel (G. moringa or G. vicinus)
% and the prey type. The difference in stomach contents between the species
% is not significant (randomization test with 100,000 replicates, P=0.11).
% 
% There are a lot of small numbers in this data set. If you pool the data
% into fish (the first five species), crustaceans (crabs and lobster), and
% octopus+unidentified
% 
% 
% 
% 

% Fish             21   18
% Crustaceans      51   14
% Misc              7    3

Z = [21 18; 51 14; 7 3];
f_chisqIndep(Z,1000)

% ans = 
% 
%       stat: 6.9443
%         df: 2
%     p_perm: 0.293
%     p_para: 0.03105

