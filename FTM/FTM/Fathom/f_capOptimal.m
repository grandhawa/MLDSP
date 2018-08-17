function result = f_capOptimal(Y,dis,X,sm,verb)
% - get optimal value of m for f_cap
%
% USAGE: result = f_capOptimal(Y,'dis',X,sm,verb);
%
% Y   = matrix of response variables (rows = obs, cols = variables)
% dis = dissimilarity measure to apply to Y
%         (e.g., dis = 'bc'; see help for f_dis)
%
% X   = (1) vector of integers specifying group membership for objects in yDis,
%       (2) ANOVA design matrix specified by dummy coding, or
%       (3) matrix of explanatory variables (rows = observations, cols = variables)
%
% sm   = use spatial median instead of centroid     (default = 0)
% verb = send output to display                     (default = 1)
%
% result = structure of results with the following fields:
%  .m     = # PCoA axes
%  .propG = corresponding proportion of yDis explained
%
%  SEE ALSO: f_cap

% -----Check inputs & set defaults:-----
if (nargin < 4), sm   = 0; end % default don't use spatial median
if (nargin < 5), verb = 1; end % default sent output to display

% -----Author:-----
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2011: propG now returns the correct value

% Get range of values of m to try:
result = f_cap(Y,dis,X,[],sm,0,0,-1,0);
m      = result.m;
clear result;

if (verb>0)
   fprintf('\nEvaluating values of m from %d to %d...\n',m(1),m(end));
end

% Evaluate m according to loo_err & loo_RSS:
nM      = numel(m);
propG   = nan(nM,1); % preallocate
loo_err = nan(nM,1);
loo_RSS = nan(nM,1);
AIC     = nan(nM,1);
evals   = zeros(nM,3);
for i = 1:nM
   cap        = f_cap(Y,dis,X,[],sm,0,0,m(i),1,[],1);
   propG(i)   = cap.Qexpl(m(i),2);
   loo_err(i) = (1-cap.loo_err.tot)*100;
   loo_RSS(i) = cap.loo_RSS;
   AIC(i)     = cap.AIC;
   evals(i,1:numel(cap.evals)) = cap.evals;
   ctr        = cap.ctr;
   method     = cap.method;
   clear cap;
end

% ----Get optimal m:---
if isequal('CDA',method)
   opt = find(loo_err==max(loo_err));
else
   opt = find(loo_RSS==min(loo_RSS));
end
opt = opt(1); % in case of ties, keep one

% % Akaike weights:
% delta = AIC - min(AIC);  % AIC differences between models
% wt    = exp(-0.5*delta); % relative likelihood of each model
% wt    = wt./nansum(wt);  % normalize, so sum to 1
% % - note: nansum needed when classifcation error = 0 and resulting AIC =
% % NaN;
%
% % Determine the best value of m:
% opt = find(wt == max(wt));            % largest wt wins
% if (length(opt)>1), opt = opt(1); end % in case of ties, select one


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   if isequal('CDA',method)
      fprintf(' Diagnostics for CAP - db-CDA:\n');
   else
      fprintf(' Diagnostics for CAP - db-CCorA: \n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf( 'm:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:\n');
   %fprintf('m  00.0000  0.0000 0.0000 0.0000 0.0000  000.00%')
   for i=1:nM
      fprintf('%d  %6.4f  %6.4f %6.4f %6.4f %6.4f  %6.2f%% \n', m(i), propG(i),...
         loo_RSS(i), evals(i,1), evals(i,2), evals(i,3), loo_err(i));
   end
   fprintf('--------------------------------------------------\n');
   fprintf('m       = # PCoA axes retained\n');
   fprintf('propG   = proportion of yDis explained by m PCoA axes\n')
   fprintf('RSS     = leave-one-out residual sums-of-squares\n')
   fprintf('d_1^2   = squared canonical correlation for axis 1\n')
   fprintf('Correct = leave-one-out classification success\n')
   fprintf('\n')
   fprintf('Central Tendency = %s \n',ctr);
   fprintf('Optimal value of m for %s may be %d\n',method, m(opt))
   fprintf('==================================================\n');
end

% -----Wrap results up into a structure:-----
result.opt   = opt; % index to optimum m
result.m     = m;
result.propG = propG;
result.AIC   = AIC;








