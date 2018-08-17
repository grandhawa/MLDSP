function result = f_ncapOptimal(Y,dis,X,grad,type,m,verb,plt)
% - evaluate a range of values of m for f_ncap
%
% USAGE: result = f_capOptimal(Y,'dis',X,'grad','type',m,verb,plt);
%
% Y    = matrix of response variables (rows = obs, cols = variables)
% dis  = dissimilarity measure to apply to Y
%        (e.g., dis = 'bc'; see help for f_dis)
%
% X    = column vector defining the univariate explanatory variable (= gradient)
%
% grad = 'von': von Bertalannfy                                (default = 'von')
%        'hyp': hyperbolic
%        'log': logistic
%        'lin': linear
%
% type = 'R2' (= default) or 'RDA'
% m    = maximum # of axes of Q to retain; m = 0 retains all       (default = 0)
% verb = send output to display                                    (default = 1)
% plt  = create diagnostic plots                                   (default = 1)
%
% result = structure of outputs with the following fields:
%  .m    = values of m evaluated
%  .stat = corresponding R2 or RDA statistic
%  .type = 'R2' or 'RDA'
%  .mOpt = optimal value if 'R2' used (= NaN if 'RDA' used)
%
%  SEE ALSO: f_ncap, f_AIC

% -----References:-----
% Anderson, D. R., K. P. Burnham, and W. L. Thompson. 2000. Null hypothesis
%  testing: problems, prevalence, and an alternative. J. Wildl. Manage.
%  64(4): 912-923.
% Burnham, K. P. and D. R. Anderson. 2001. Kullback-Leibler information as
%  a basis for strong inference in ecological studies. Wildlife Research 28:
%  111-119.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. pages 612-616
% Millar, R. B., M. J. Anderson, & G. Zunun. 2005. Fitting nonlinear
%  environmental gradients to community data: a general distance-based approach.
%  Ecology 86(8): 2245-2251.
%
% Portions of this function are based on ideas implemented in Millar's
% 'NCAP.R' code for R.

% -----Author:-----
% by David L. Jones, Mar-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check inputs & set defaults:-----
if (nargin < 4), grad = 'von'; end % default send results to display
if (nargin < 5), type = 'R2';  end % default stat is R^2
if (nargin < 6), m    = 0';    end % default use all axes of Q
if (nargin < 7), verb = 1;     end % default sent output to display
if (nargin < 8), plt  = 1;     end % default create plots
% --------------------------------------

% Get range of values of m to try:
opt = 1; % indicate f_ncap is called by f_ncapOptimum
if (m==0)
   result = f_ncap(Y,dis,X,grad,type,0,(0),0,0,opt);
else
   result = f_ncap(Y,dis,X,grad,type,0,(m),0,0,opt);
end
m      = (1:result.m)'; % range of m's to evaluate
propY  = result.propY;  % cumulative % of Y explained by axes
clear result;

if (verb>0)
   fprintf('\nEvaluating values of m from %d to %d...\n',m(1),m(end));
end

switch upper(type)
   case 'R2'
      % Evaluate m according to an R^2-based AIC:
      nM    = numel(m);
      b     = nan(nM,1); % preallocate
      stat  = nan(nM,1);
      AIC   = nan(nM,1);
      for i = 1:numel(m)
         ncap    = f_ncap(Y,dis,X,grad,'R2',0,m(i),0,0,0);
         b(i)    = ncap.b;
         stat(i) = ncap.stat;
         AIC(i)  = ncap.AIC;
         clear ncap;
      end
      
      % Akaike weights:
      delta = AIC - min(AIC);  % AIC differences between models
      wt    = exp(-0.5*delta); % relative likelihood of each model
      wt    = wt./sum(wt);     % normalize, so sum to 1
      
      % Determine the best value of m:
      idx = find(wt == max(wt));            % largest wt wins
      if (length(idx)>1), idx = idx(1); end % in case of ties, select one
      mOpt = m(idx);
      
      % Show optimal m:
      if (verb>0)
         fprintf('\nThe optimal value of m (between %d & %d) = %d\n',m(1),m(end),mOpt);
      end
      
   case 'RDA'
      % Evaluate m according to the RDA statistic:
      nM    = numel(m);
      b     = nan(nM,1); % preallocate
      stat  = nan(nM,1);
      for i = 1:numel(m)
         ncap    = f_ncap(Y,dis,X,grad,'RDA',0,m(i),0,0,0);
         b(i)    = ncap.b;
         stat(i) = ncap.stat;
         clear ncap;
      end
      mOpt = NaN;
   otherwise
      error('TYPE should be ''R2'' or ''RDA''!');
end


% Optionally plot results:
if (plt>0)
   % -----Plot M vs. Percent Variation Explained:-----
   figure('Name','NCAP: Percent Variation Explained');
   set(gcf,'color','w'); % set background color to white
   plot(1:numel(propY),propY,'b-');
   xlabel('m');
   ylabel('% Variation Explained');
   set(gca,'Xgrid','on');
   
   % -----Plot M vs. b-hat:-----
   figure('Name','NCAP: Regression Coefficients');
   set(gcf,'color','w'); % set background color to white
   hold on;
   plot(m,b,'bo');       % regression coefficients
   if isequal(type,'R2'), plot(mOpt,b(idx),'r*'); end
   xlabel('m');
   ylabel('b-hat');
   box on;
   set(gca,'Xlim',[0 m(end)*1.05]);
   
   switch type
      case 'R2'
         % -----Plot M vs. R2 Statistic:-----
         figure('Name','NCAP: Optimal m');
         set(gcf,'color','w'); % set background color to white
         hold on;
         plot(m,stat,'bo');      % observed R2 stat
         plot(mOpt,stat(idx),'r*');
         txt = sprintf('Optimal m = %d',mOpt);
         title(txt);
         xlabel('m');
         ylabel('Correlation (R^2)');
         box on;
         set(gca,'Xlim',[0 m(end)*1.05],'Ylim',[0 1]);
         
      case 'RDA'
         % -----Plot M vs. RDA Statistic:-----
         figure('Name','NCAP: m vs. RDA statistic');
         set(gcf,'color','w'); % set background color to white
         hold on;
         plot(m,stat,'bo');      % observed RDA stat
         xlabel('m');
         ylabel('RDA statistic');
         box on;
         set(gca,'Xlim',[0 m(end)*1.05]);
   end
end


% Wrap results up into a structure:
result.m    = m;    % values of m evaluated
result.stat = stat; % corresponding R2 or RDA statistic
result.type = type; % 'R2' or 'RDA'
result.mOpt = mOpt; % optimal value if 'R2' used

