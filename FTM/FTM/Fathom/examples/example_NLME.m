% Example of fitting a Nonlinear Mixed Effects model
% by David L. Jones, Aug-2014

% The file 'orange_trees.csv' contains data concerning the growth of orange
% trees from Draper and Smith (1981) as follows:
% 
% grp = individual tree
% cir = circumference (mm)
% day = the elapsed day each measurement was taken


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    NOTES:                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This example includes steps for calculating F-ratio's and permuted p-values
% to test the statistical significance of the output of the 'nlmefitsa'
% function. This is used instead of the regular 'nlmefit' function as the latter
% does not work when the data are permuted. You will need to increase the value
% of 'iter' from 10 to, say, 1000 to get a proper p-value, but it may take some
% time to run as the fitting procedure used by 'nlmefitsa' is pretty computer
% intensive.
% 
% As an alternative to using 'nlmefitsa' to calculate a permuted p-value to
% assess the significance of the corresponding F-ratio, you could try using the
% regular 'nlmefit' function and compare its AIC value against one produced from
% a null model, say, a regression that produces just a horizontal line. In that
% case, you could use the Information-Theoretic approach to say the nonlinear
% regression is 'substantially' better than a null model (vs. saying it's
% 'significantly' better when assessing an F-ratio via a p-value). You might get
% some ideas concerning the 'horizontal line' model here:      
% 
% http://www.mathworks.com/matlabcentral/answers/27484-p-value-from-nlinfit

% Clear the workspace:
clz;

% Load the data:
raw = f_importCSV('orange_trees.csv',1)

% Parse the data:
grp = raw.dat(:,1); % group (= individual tree)
day = raw.dat(:,2); % day
cir = raw.dat(:,3); % circumference
clear raw;

% Plot the data by group:
uGrp = f_unique(grp); % get list of unique groups
nGrp = numel(uGrp);   % get # unique groups
figure; 
set(gcf,'color','w');
hold on;
for i = 1:nGrp
   idx = grp==uGrp(i); % get index to rows for this group
   h(i) = plot(day(idx),cir(idx),'o','Color',f_rgb(i),'LineWidth',2);
end
xlabel('Time (days)')
ylabel('Circumference (mm)')
title('{\bf Orange Tree Growth}')
legend([repmat('Tree ',nGrp,1),num2str(uGrp)],'Location','NW')
grid on;
box on;

% Specify a 3 parameter logistic growth function:
fun = @(PHI,t)(PHI(:,1))./(1+exp(-(t-PHI(:,2))./PHI(:,3)));

% Define initial model parameters:
beta0 = [100 100 100]';

% Fit a nonlinear mixed effects model:
[beta,PSI,stats,B] = nlmefitsa(day,cir,grp,[],fun,beta0);

beta
PSI
stats
B

% beta =
% 
%        194.29
%        739.94
%         354.9
% 
% 
% PSI =
% 
%        1024.9            0            0
%             0       47.522            0
%             0            0       71.536
% 
% 
% stats = 
% 
%           logl: []
%            aic: []
%            bic: []
%         sebeta: []
%            dfe: 28
%           covb: []
%     errorparam: 7.8624
%           rmse: 8.0541
%           ires: [35x1 double]
%           pres: [35x1 double]
%          iwres: [35x1 double]
%          pwres: [35x1 double]
%          cwres: [35x1 double]
% 
% 
% B =
% 
%       -29.937       31.816      -37.653       40.254      -5.3488
%       -1.2828     -0.63408     -0.38941      0.23149       1.4557
%        2.6647      -0.1398       1.5667      -2.5329      -1.5856


% Use the estimated fixed effects in 'beta' and the estimated random effects for
% each tree in 'B' to plot the model through the data:
PHI   = repmat(beta,1,nGrp) + B; % Fixed effects + Random effects
tplot = 0:0.1:1600;
for i = 1:nGrp
  fitted_model = @(t)(PHI(1,i))./(1+exp(-(t-PHI(2,i))./ PHI(3,i)));
  plot(tplot,fitted_model(tplot),'Color',f_rgb(i),'LineWidth',2)
end

% Calculate ANOVA statistics for the regression:
Y    = cir;                % define response variable
Yres = stats.ires;         % residuals
SSt  = trace(Y'*Y);        % Sum-of-Squares total
SSe  = trace(Yres'*Yres);  % Sum-of-Squares error
SSr  = SSt - SSe;          % Sum-of-Squares regression
dfE  = stats.dfe;          % Degrees-of-Freedom error (= N - # model parameters)
dfR  = size(Yres,1) - dfE; % Degrees-of-Freedom regression (= # model parameters)
MSe  = SSe/dfE;            % Mean Squares error (= stats.rmse^2)
MSr  = SSr/dfR;            % Mean Squares regression model 
F    = MSr/MSe             % F-ratio
R2   = 1 - SSe/SSt         % R-squared

% F  = 1262.4
% R2 = 0.99684

% -----Permutation Test:-----
iter    = 10;
F_perm  = zeros(iter-1,1); % preallocate result array
fprintf('\nCalculating permuted F-ratios:\n');
for i = 1:(iter-1) % observed value is considered a permutation
   % Display progress:
   fprintf('%d of %d...\n',i,iter-1);
   
   % Close figure windows:
   close all;
   
   [~,~,permStats] = nlmefitsa(f_shuffle(day),cir,grp,[],fun,beta0);
   
   Yres_perm = permStats.ires;               % residuals
   SSe_perm  = trace(Yres_perm'*Yres_perm);  % Sum-of-Squares error
   SSr_perm  = SSt - SSe_perm;               % Sum-of-Squares regression
   dfE_perm  = permStats.dfe;                % Degrees-of-Freedom error (= N - # model parameters)
   dfR_perm  = size(Yres_perm,1) - dfE_perm; % Degrees-of-Freedom regression (= # model parameters)
   MSe_perm  = SSe_perm/dfE_perm;            % Mean Squares error (= stats.rmse^2)
   MSr_perm  = SSr_perm/dfR_perm;            % Mean Squares regression model
   F_perm(i) = MSr_perm/MSe_perm;            % permuted F-ratio
end
p = (sum(F_perm>=F)+1)./(iter); % count values & convert to probability
% ---------------------------
























