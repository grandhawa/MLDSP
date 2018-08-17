% Convert sex and fat to one term (interaction term) to perform the 'Preliminary
% ANOVA' in Sokal & Rohlf, p. 326:


[null,dave] = f_designMatrix([sex fat])

% null =
% 
%      1     0     1
%      1     0     1
%      1     0     1
%      1     0     0
%      1     0     0
%      1     0     0
%      0     1     1
%      0     1     1
%      0     1     1
%      0     1     0
%      0     1     0
%      0     1     0
% 
% 
% dave =
% 
%      1     0     0
%      1     0     0
%      1     0     0
%      0     1     0
%      0     1     0
%      0     1     0
%      0     0     1
%      0     0     1
%      0     0     1
%      0     0     0
%      0     0     0
%      0     0     0


grp = [1 1 1 2 2 2 3 3 3 4 4 4]'


% grp =
% 
%      1
%      1
%      1
%      2
%      2
%      2
%      3
%      3
%      3
%      4
%      4
%      4
% 
% dum = f_designMatrix(grp)
% 
% dum =
% 
%      1     0     0
%      1     0     0
%      1     0     0
%      0     1     0
%      0     1     0
%      0     1     0
%      0     0     1
%      0     0     1
%      0     0     1
%      0     0     0
%      0     0     0
%      0     0     0
% 
% f_npManova(dis,[grp],1000,1);
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'       'MS'        'F'         'p'    
%     'factor 1'    [ 3]    [65904]    [ 21968]    [15.064]    [0.003]
%     'residual'    [ 8]    [11667]    [1458.3]    [   NaN]    [  NaN]
%     'total'       [11]    [77570]    [   NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------


