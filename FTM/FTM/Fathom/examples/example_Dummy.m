% Example of dummy coding categorical variables for regression analysis
% 
% by David L. Jones, Oct-2010
% 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Create some categorical variables:
var_01 = [ones(12,1);ones(12,1)+1]              % 1 factor coding for 2 levels (n = 24)
var_02 = repmat([1 1 1 1 2 2 2 2 3 3 3 3]',2,1) % 1 factor coding for 3 levels (n = 24)

% -----Create dummy codes using f_dummy:-----
help f_dummy
% - dummy coding of categorical variables
%  
%   USAGE: X = f_dummy(Y,trim)
%  
%   Y    = column vector of integers specifying factor levels or group membership;
%   trim = trim last column to avoid a singular matrix (default = 1)
%   
%   SEE ALSO: f_dummy2cat, f_xMatrix, f_helmert, f_designMatrix, f_modelMatrix

var_01_Rx = f_dummy(var_01)
% var_01_Rx =
% 
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0

var_02_Rx = f_dummy(var_02)
% var_02_Rx =
% 
%      1     0
%      1     0
%      1     0
%      1     0
%      0     1
%      0     1
%      0     1
%      0     1
%      0     0
%      0     0
%      0     0
%      0     0
%      1     0
%      1     0
%      1     0
%      1     0
%      0     1
%      0     1
%      0     1
%      0     1
%      0     0
%      0     0
%      0     0
%      0     0


% -----Create dummy codes using f_xMatrix:-----
% help f_xMatrix
%   - design matrix of contrast codes for ANOVA linear models
%  
%   USAGE: [X,int] = f_xMatrix(Y,orth,nest)
%  
%   Y    = matrix of integers specifying factor levels or group membership;
%          (rows = observations, cols = factors)
%   orth = create orthogonal contrast codes                       (default = 0)
%   nest = factors are nested (i.e., factor 2 within factor 1)    (default = 0)
%          Eaxmple:
%          factor 1 = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]'; % main factor
%          factor 2 = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]'; % nested factor
%  
%   X    =  design matrix (all crossed factors combined, or 1 nested factor)
%   int  = interaction term(s) for un-nested factors
%   
%   SEE ALSO: f_dummy, f_helmert, f_designMatrix, f_modelMatrix

var_01_Rx = f_xMatrix(var_01)
% var_01_Rx =
% 
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%      1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1


var_02_Rx = f_xMatrix(var_02)
% var_02_Rx =
% 
%      1     0
%      1     0
%      1     0
%      1     0
%     -1     1
%     -1     1
%     -1     1
%     -1     1
%      0    -1
%      0    -1
%      0    -1
%      0    -1
%      1     0
%      1     0
%      1     0
%      1     0
%     -1     1
%     -1     1
%     -1     1
%     -1     1
%      0    -1
%      0    -1
%      0    -1
%      0    -1


% -----Create the 'var_01 x var_02' interaction term:-----
[null,var_01_02_Rx] = f_xMatrix([var_01 var_02]);
var_01_02_Rx
% var_01_02_Rx =
% 
%      1     0
%      1     0
%      1     0
%      1     0
%     -1     1
%     -1     1
%     -1     1
%     -1     1
%      0    -1
%      0    -1
%      0    -1
%      0    -1
%     -1     0
%     -1     0
%     -1     0
%     -1     0
%      1    -1
%      1    -1
%      1    -1
%      1    -1
%      0     1
%      0     1
%      0     1
%      0     1

