% Examples of Canonical Analysis of Principal Coordinates
% 
% by David L. Jones, Feb-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% Apr-2011: updated
% Nov-2012: added Fisher's Iris data example
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Bristol Channel Plankton:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% File: '.../examples/bristol.mat'
% 
% The file 'bristol.mat' contains data from Collins & Williams (1982)
% concerning abundances (# per cubic meter) of 24 species of holozooplankton
% species collected from double oblique plankton net hauls at 57 sites in the
% Bristol Channel, UK. The 57 sites were grouped according to the following 4
% categories describing the salinity regime: TEst = true estuarine, EstM =
% estuarine and marine, Sten = stenohaline marine, and Eury = euryhaline marine.
% Note there is no site 30.

% Load the data:
load bristol.mat;

% 4th-root transform data:
zoo_4 = f_normal(zoo,'4');

% Find optimal value of m:
f_capOptimal(zoo_4,'bc',grp,0,1);
% 
% Evaluating values of m from 2 to 7...
% ==================================================
%  Diagnostics for CAP - db-CDA:
% --------------------------------------------------
% m:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:
% 2  69.5416  3.9676 0.9264 0.7549 0.0000   94.74% 
% 3  78.5365  3.8044 0.9293 0.7866 0.5517   96.49% 
% 4  85.0543  3.8320 0.9327 0.7998 0.5517   96.49% 
% 5  89.5931  3.9399 0.9327 0.8358 0.5905   98.25% 
% 6  93.5448  3.9099 0.9332 0.8399 0.6643  100.00% 
% 7  97.1428  4.0540 0.9380 0.8498 0.6685  100.00% 
% --------------------------------------------------
% m       = # PCoA axes retained
% propG   = proportion of yDis explained by m PCoA axes
% RSS     = leave-one-out residual sums-of-squares
% d_1^2   = squared canonical correlation for axis 1
% Correct = leave-one-out classification success
% 
% Central Tendency = centroid 
% Optimal value of m for CDA may be 6
% ==================================================


% Perform CAP using optimal m:
cap = f_cap(zoo_4,'bc',grp,[],0,1000,1,6,1);
% 
% ==================================================
%  CAP - Canonical Discriminant Analysis:
% --------------------------------------------------
% Trace Stat    = 2.4373  p =  0.00100 
% Greatest Root = 0.9332  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% No. of axes of Q used (m)     = 6 
% Variability of yDis explained = 93.54 % 
% Canonical Correlations:
%   0.9660  0.9165  0.8150
% Squared Canonical Correlations (= delta^2):
%   0.9332  0.8399  0.6643
% ==================================================
% 
% ==================================================
%             LOO CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Correct  
%    1            100.0 % 
%    2            100.0 % 
%    3            100.0 % 
%    4            100.0 % 
% 
% 
% Total Correct  = 100.00 % 
% Total Error    = 0.00 % 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      3      4      
%      1  100.0    0.0    0.0    0.0 
%      2    0.0  100.0    0.0    0.0 
%      3    0.0    0.0  100.0    0.0 
%      4    0.0    0.0    0.0  100.0 
% 
% ==================================================

% Test the significance of the observed classification success rate:
f_chanceClass(grp,1,1000,1);
%
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1             3.08 % 
%    2             9.97 % 
%    3             6.93 % 
%    4             6.03 % 
% 
% Total Correct = 26.01 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 25.97 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% -> Reject the null hypothesis that the observed re-classification success
%    rate is no better than that expected by chance alone, at alpha = 0.05.   

% Create plot:
f_capPlot(cap,f_unique(grp_txt),[],zoo_4,zoo_txt_small,0.05,'none',1,1);

% Save graphics as PDF's:
f_pdf('bristol_db-CDA')
f_pdf('bristol_speciesImp')
f_pdf('bristol_WAscores')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Fisher's Iris data:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% File: '.../examples/iris.mat'
% 
% This example uses a subset of Fisher's data to demonstrate plotting of
% jittered axes when only a single canonical axis is present.

% Load the data:
load iris.mat

% Create index to just groups 1 and 2:
idx = find(grps>=1 & grps<=2);

% Find optimal value of m:
f_capOptimal(iris(idx,:),'euc',grps(idx),1,1);
% 
% Evaluating values of m from 1 to 3...
% 
% ==================================================
%  Diagnostics for CAP - db-CDA:
% --------------------------------------------------
% m:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:
% 1  90.5393  0.0657 0.9356 0.0000 0.0000  100.00% 
% 2  97.9848  0.0408 0.9615 0.0000 0.0000  100.00% 
% 3  99.6582  0.0412 0.9619 0.0000 0.0000  100.00% 
% --------------------------------------------------
% m       = # PCoA axes retained
% propG   = proportion of yDis explained by m PCoA axes
% RSS     = leave-one-out residual sums-of-squares
% d_1^2   = squared canonical correlation for axis 1
% Correct = leave-one-out classification success
% 
% Central Tendency = spatial median 
% Optimal value of m for CDA may be 1
% ==================================================

% Perform CAP using optimal m:
cap = f_cap(iris(idx,:),'euc',grps(idx),[],1,1000,1,1,1);
% 
% ==================================================
%  CAP - Canonical Discriminant Analysis:
% --------------------------------------------------
% Trace Stat    = 0.9356  p =  0.00100 
% Greatest Root = 0.9356  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% No. of axes of Q used (m)     = 1 
% Variability of yDis explained = 90.54 % 
% Canonical Correlations:
%   0.9672
% Squared Canonical Correlations (= delta^2):
%   0.9356
% ==================================================
% 
% ==================================================
%             LOO CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Correct  
%    1            100.0 % 
%    2            100.0 % 
% 
% 
% Total Correct  = 100.00 % 
% Total Error    = 0.00 % 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      
%      1  100.0    0.0 
%      2    0.0  100.0 
% 
% ==================================================

% Test the significance of the observed classification success rate:
f_chanceClass(grps(idx),1,1000,1);
% ==================================================
%            PROPORTIONAL CHANCE CRITERION:
% --------------------------------------------------
% Group        Correct  
%    1            25.00 % 
%    2            25.00 % 
% 
% Total Correct = 50.00 % 
% --------------------------------------------------
% 
% Mean randomized classification success = 50.36 % 
% p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% -> Reject the null hypothesis that the observed re-classification success
%    rate is no better than that expected by chance alone, at alpha = 0.05.
    
% Create plot:
 f_capPlot(cap,{'grp 1' 'grp 2'},[],iris(idx,:),[],0.05,'none',1,1);

 % Save graphics as PDF's:
f_pdf('iris_db-CDA')
f_pdf('iris_Imp')
f_pdf('iris_WAscores')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Classify Unknowns:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% File: '.../examples/iris.mat'

% Load the data:
load iris.mat

% Randomly allocate the data into baseline (= TRAINING) and mixture (= UNKNOWN) sets:
idx  = f_shuffle([1:numel(grps)]');
idxB = idx(1:140);   % index to baseline samples
idxM = idx(141:end); % index to mixture samples
B    = iris(idxB,:); % baseline data
Bgrp = grps(idxB);   % group membership of baseline data
M    = iris(idxM,:); % mixture data
Mgrp = grps(idxM);   % known group membership of mixture data

% Get index to sort rows by group:
[nul,idxS] = sort(Bgrp);

% Find optimal value of m:
f_capOptimal(B(idxS,:),'euc',Bgrp(idxS),1,1); % <- **it's critical data are sorted by group!**
% 
% Evaluating values of m from 1 to 4...
% 
% ==================================================
%  Diagnostics for CAP - db-CDA:
% --------------------------------------------------
% m:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:
% 1  92.2517  1.4748 0.9283 0.0000 0.0000   91.43% 
% 2  97.6777  1.4162 0.9637 0.0943 0.0000   95.00% 
% 3  99.5034  1.5302 0.9672 0.2099 0.0000   93.57% 
% 4  100.0000  1.5410 0.9680 0.2169 0.0000   94.29% 
% --------------------------------------------------
% m       = # PCoA axes retained
% propG   = proportion of yDis explained by m PCoA axes
% RSS     = leave-one-out residual sums-of-squares
% d_1^2   = squared canonical correlation for axis 1
% Correct = leave-one-out classification success
% 
% Central Tendency = spatial median 
% Optimal value of m for CDA may be 2
% ==================================================


% Perform CAP using optimal m:
cap = f_cap(B(idxS,:),'euc',Bgrp(idxS),[],1,1000,1,2,1);
% 
% ==================================================
%  CAP - Canonical Discriminant Analysis:
% --------------------------------------------------
% Trace Stat    = 1.0580  p =  0.00100 
% Greatest Root = 0.9637  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% No. of axes of Q used (m)     = 2 
% Variability of yDis explained = 97.68 % 
% Canonical Correlations:
%   0.9817  0.3070
% Squared Canonical Correlations (= delta^2):
%   0.9637  0.0943
% ==================================================
% 
% ==================================================
%             LOO CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Correct  
%    1            100.0 % 
%    2             89.4 % 
%    3             95.7 % 
% 
% 
% Total Correct  = 95.00 % 
% Total Error    = 5.00 % 
% 
% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      3      
%      1  100.0    0.0    0.0 
%      2    0.0   89.4   10.6 
%      3    0.0    4.3   95.7 
% 
% ==================================================

% Classify observations in mixture (UNKNOWN) set:
capM = f_cap(B(idxS,:),'euc',Bgrp(idxS),M,1,0,0,2,0); 

% Compare predicted group membership with actual group membership:
err = f_errRate(Mgrp,capM.grpW);
fprintf('Total Correct  = %4.2f %s \n',(1-err.tot)*100,'%');
% Total Correct  = 100.00 %

