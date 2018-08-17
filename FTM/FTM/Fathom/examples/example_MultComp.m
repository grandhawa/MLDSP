% Examples of Multiple Comparison Tests
% 
% by David L. Jones, Oct-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/glass.mat'
% These data consist of USA Forensic Science Service measurements of the
% oxide content from 214 samples of 6 types of glass. These data were
% obtained from the UCI Machine Learning:
% http://archive.ics.uci.edu/ml/datasets/Glass+Identification

% Variables:
% Y     = matrix of response variables (either refractive index or % weight of 8
%         elements in their corresponding oxides)
% Y_txt = cell array of corresponding labels (= RI, Na, Mg, Al, Si, K, Ca, Ba, Fe) 
% grp   = column vector of integers specifying type of glass:
%         1: building_windows_float_processed
%         2: building_windows_non_float_processed
%         3: vehicle_windows_float_processed
%         4: containers
%         5: tableware
%         6: headlamps

% Load the data:
load glass.mat;

% Create Euclidean distance matrix:
yDis = f_dis(Y,'euc');

% Test the null hypothesis of no significant difference in refractive index
% and oxide levels among six types of glass:
f_npManova(yDis,grp,1000,1);
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'     'SS'        'MS'        'F'         'p'    
%     'factor 1'    [  5]    [431.55]    [86.311]    [19.702]    [0.001]
%     'residual'    [208]    [ 911.2]    [4.3808]    [   NaN]    [  NaN]
%     'total'       [213]    [1342.8]    [   NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% (Note: NaNs are placeholders for the ANOVA table)
% 
% (Data with replication)
% 
% -> reject the null hypothesis at an alpha = 0.05 level.

% Calculate multivariate dispersions:
md = f_npDisp(yDis,grp,1,1);
% ==================================================
%       Quantify Multivariate Dispersions:
% ==================================================
% 
% # Pos Eigenvalues = 9
% # Neg Eigenvalues = 0
% 
% Average distance to spatial median:
% Group 1 = 0.9266
% Group 2 = 1.6531
% Group 3 = 0.8106
% Group 4 = 2.8999
% Group 5 = 1.9847
% Group 6 = 1.6507
% ---------------------------------------------------

% Test for significant differences in among-group dispersion:
f_npManova(f_dis(md.z,'euc'),grp,1000,1);


% >>>>> Might need to transform these, or standardize, or drop the RI to
% homogenize the variances.






