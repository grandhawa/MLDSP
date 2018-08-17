% Example of calculating bootstrapped confidence intervals
% by David L. Jones, 12-Sep-2013

% Create some data:
X = repmat([1:25]',1,10)
% X =
% 
%      1     1     1     1     1     1     1     1     1     1
%      2     2     2     2     2     2     2     2     2     2
%      3     3     3     3     3     3     3     3     3     3
%      4     4     4     4     4     4     4     4     4     4
%      5     5     5     5     5     5     5     5     5     5
%      6     6     6     6     6     6     6     6     6     6
%      7     7     7     7     7     7     7     7     7     7
%      8     8     8     8     8     8     8     8     8     8
%      9     9     9     9     9     9     9     9     9     9
%     10    10    10    10    10    10    10    10    10    10
%     11    11    11    11    11    11    11    11    11    11
%     12    12    12    12    12    12    12    12    12    12
%     13    13    13    13    13    13    13    13    13    13
%     14    14    14    14    14    14    14    14    14    14
%     15    15    15    15    15    15    15    15    15    15
%     16    16    16    16    16    16    16    16    16    16
%     17    17    17    17    17    17    17    17    17    17
%     18    18    18    18    18    18    18    18    18    18
%     19    19    19    19    19    19    19    19    19    19
%     20    20    20    20    20    20    20    20    20    20
%     21    21    21    21    21    21    21    21    21    21
%     22    22    22    22    22    22    22    22    22    22
%     23    23    23    23    23    23    23    23    23    23
%     24    24    24    24    24    24    24    24    24    24
%     25    25    25    25    25    25    25    25    25    25

% Calculate bootstrapped 95% confidence intervals from 1000 resamples:
result = f_bootCI(X,0.95,1000)
% result = 
% 
%     mean: [13 13 13 13 13 13 13 13 13 13]
%     ci_u: [15.4 15.44 15.48 15.38 15.46 15.42 15.3 15.34 15.2 15.28]
%     ci_l: [10.6 10.56 10.52 10.62 10.54 10.58 10.7 10.66 10.8 10.72]
%     conf: 0.95

% Use Matlab's 'bootci' function:
myfun = @(x)mean(x);
mlab  = bootci(1000,{myfun,X})
% mlab =
% 
%   Columns 1 through 7
% 
%         10.24        10.24        10.24        10.24        10.24        10.24        10.24
%         15.68        15.68        15.68        15.68        15.68        15.68        15.68
% 
%   Columns 8 through 10
% 
%         10.24        10.24        10.24
%         15.68        15.68        15.68    
% 
% -> unlike 'f_bootCI', Matlab's 'bootci' function does not create bootstrapped
%    samples SEPARATELY for each column of the input data.
