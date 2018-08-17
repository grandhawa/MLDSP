% Example of calculating Species Indicator Values
% 
% by David L. Jones, Apr-2013
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/'Dufrene_Legendre.csv'
% File: '.../examples/'Legendre_Legendre_p400.csv'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Example: Dufrene & Legendre, 1997                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% File: 'Dufrene_Legendre.csv' is example data from Table 2 of Dufrene &
% Legendre, 1997.

% Load data:
raw = f_importCSV('Dufrene_Legendre.csv');

% Parse the data:
Y     = raw.dat(:,3:5);
grp   = raw.dat(:,1);
site  = raw.dat(:,2);
Y_txt = {'sp_A' 'sp_B' 'sp_C'};
clear raw;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Cluster of 1 group:          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Create a 1-cluster grouping variable:
grp1 = ones(25,1);
% 
% Calculate indicator values:
f_indVal(Y,grp1,0)
% ans = 
% 
%     indVal: [100 72 40]
%          p: NaN
%         IV: [100 72 40]
%          S: [1 1 1]
%          G: [1 1 1]
% 
% -> Note that since there is only one group here there is no corresponding
%    statistical test, so the number of permutation iterations in the
%    f_indVal function must be set to 0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Cluster of 2 groups:          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Create a 2-cluster grouping variable:
grp2 = [ones(15,1); ones(10,1)*2];
% 
% Calculate indicator values:
f_indVal(Y,grp2,1000)
% ans = 
% 
%     indVal: [60.87 85.714 66.667]
%          p: [0.001 0.001 0.003]
%         IV: [2x3 double]
%          S: [2x3 logical]
%          G: [1 1 1]
% 
% -> All indicator values (indVal) are statistically significant at the
%    alpha = 0.05 level.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Example: Legendre & Legendre, 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% File: 'Legendre_Legendre.csv' is example data from Table 8.10 of Legendre &
% Legendre, 2012.

% Clear the workspace:
clz;

% Load data:
raw = f_importCSV('Legendre_Legendre_p400.csv');

% Parse the data:
Y     = raw.dat(:,3:5);
grp   = raw.dat(:,1);
site  = raw.dat(:,2);
Y_txt = {'sp_A' 'sp_B' 'sp_C'};
clear raw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Cluster of 5 groups:          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Calculate indicator values:
f_indVal(Y,grp,1000)
% 
% ans = 
% 
%     indVal: [30 40 90]
%          p: [0.001 0.001 0.001]
%         IV: [5x3 double]
%          S: [5x3 logical]
%          G: [3 1 1]
% 
% -> All indicator values (indVal) are statistically significant at the
%    alpha = 0.05 level.
% 
% Note: the indVal.G variable indicates that species 1 is an indicator for
%       group 3 and the other species are indicators for group 1. Species 3
%       has the the greatest indicator power at 90%.
