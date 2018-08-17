% Distributional properties of Greenwood's statistic
% by David L. Jones, Apr-2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 1:                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clz;

% Set the sample size:
n = 1;

% Show random distribution:
R = f_greenwood_rnd(n,1000*10,1);

% Plot probability densities:
[G,PD] = f_greenwood_plt(n,'PDF');
[~,CD] = f_greenwood_plt(n,'CDF');

sum(R.G==0.5)
% ans = 0 -> very unlikely to get such an extreme value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 2:                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clz;

% Set the sample size:
n = 2;

% Show random distribution:
f_greenwood_rnd(n,1000*10,1);

% Calculate probability densities:
G  = ((1/(n+1)):0.01:1)';
PD = f_greenwood_pdf(G,n);
CD = f_greenwood_cdf(G,n);

% Plot probability densities:
[G,PD] = f_greenwood_plt(n,'PDF');
[~,CD] = f_greenwood_plt(n,'CDF');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 3:                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clz;

% Set the sample size:
n = 3;

% Show random distribution:
f_greenwood_rnd(n,1000*10,1);

% Calculate probability densities:
G  = ((1/(n+1)):0.01:1)';
PD = f_greenwood_pdf(G,n);
CD = f_greenwood_cdf(G,n);

% Plot probability densities:
[G,PD] = f_greenwood_plt(n,'PDF');
[~,CD] = f_greenwood_plt(n,'CDF');

% Re-create Table 1 of Gardner (1952):
x = [0.25 0.2875 0.325 0.333 0.4 0.475 0.5 0.525 0.55 0.6 0.614 0.7 0.8 0.9 1]';
y = f_greenwood_pdf(x,n);
p = f_greenwood_cdf(x,n);
[x y p]
% ans =
%          0.25            0            0
%        0.2875       3.6502     0.091255
%         0.325       5.1622      0.25811
%         0.333       5.4305      0.30049
%           0.4       3.5824      0.60008
%         0.475       1.9417      0.80516
%           0.5        1.458       0.8476
%         0.525       1.1235      0.87954
%          0.55      0.89241      0.90459
%           0.6      0.57847      0.94074
%         0.614      0.51337      0.94837
%           0.7      0.24113      0.97951
%           0.8     0.085142      0.99483
%           0.9     0.017622      0.99944
%             1   1.7764e-15            1


% Get parameters of the distribution for n=3:
par = f_greenwood_par(3);                                % theoretical
R    = f_greenwood_rnd(n,1000*100,1);                    % random simulation
sim  = [mean(R.G) std(R.G) skewness(R.G) kurtosis(R.G)]; % Monte Carlo

% Comare theoretical vs. simulated:
[par' sim']
% ans =
% 
%           0.4      0.39918
%        0.1069      0.10609
%        1.3512       1.3396
%        5.1439       5.1024



 
% ans =
% 
%             5       1.5869       6.8274
%            10       1.7065       8.3506
%            50       1.2184       6.3779
%           100      0.92619       5.0256
%          1000      0.31375       3.2411
