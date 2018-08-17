% Example of using Greenwood's Statistic
% 
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Gene sequence data:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This example is from:
% Riley, M. C., A. Clare, and R. D. King. 2007. Location distribution of gene
% functional classes in Arabidopsis thaliana. BMC Bioinformatics 8:112
% 
% http://www.biomedcentral.com/1471-2105/8/112

% Clear workspace:
clz;

% X is a sequence 55 units long consisting of 11 evenly spaced values that
% differ by 5.5 units:
X = (0:5.5:55)';
f_greenwood(X,1,2,1000);
% 
% ==================================================
%               Greenwood's Statistic:             
% --------------------------------------------------
% G stat  = 0.1000 (n = 9) 
% p-value = 1.00000 
% p_perm  = 1.00000 (with 1000 iterations)
% --------------------------------------------------

% Y is a sequence the same length with the first 6 values evenly spaced by 10
% units, and the last 5 values spaced only 1 unit apart:
Y = [0:10:50 51:55]';
f_greenwood(Y,1,2,1000);
% 
% ==================================================
%               Greenwood's Statistic:             
% --------------------------------------------------
% G stat  = 0.1669 (n = 9) 
% p-value = 0.55677 
% p_perm  = 0.57100 (with 1000 iterations)
% --------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Bus station arrivals:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This example is from:
% http://en.wikipedia.org/wiki/Greenwood_statistic

% Clear workspace:
clz;

% A depicts the arrival times of 11 buses at a station during a 1 hour period,
% where each arrives 6 minutes apart:
A = (0:6:60)';
f_greenwood(A,1,1,1000);
% 
% ==================================================
%               Greenwood's Statistic:             
% --------------------------------------------------
% G stat  = 0.1000 (n = 9) 
% p-value = 1.00000 
% p_perm  = 1.00000 (with 1000 iterations)
% --------------------------------------------------

% In scenario B, 6 buses arrived 10 minutes apart with each of the remaining 5
% arriving 2 minutes apart.
B = [0:10:50 52:2:60]';
f_greenwood(B,1,1,1000);
% 
% ==================================================
%               Greenwood's Statistic:             
% --------------------------------------------------
% G stat  = 0.1444 (n = 9) 
% p-value = 0.83250 
% p_perm  = 0.84400 (with 1000 iterations)
% ---------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Miscellaneous Examples:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Clear the workspace:
clz;

% Create data:
Z = [19 3 3 1]';
f_greenwood(Z,1,1,1000);
% 
% ==================================================
%               Greenwood's Statistic:             
% --------------------------------------------------
% G stat  = 0.8025 (n = 2) 
% p-value = 0.03540 
% p_rand  = 0.03500 (with 1000 iterations)
% --------------------------------------------------
