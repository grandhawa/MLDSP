function [G,n,p,p_rand] = f_greenwood(X,verb,plt,iter)
% - Greenwood's statistic
%
% USAGE: [G,n,p,p_rand] = f_greenwood(X,verb,plt,iter)
%
% X    = column vector of input data
% verb = optionaly send results to display                         (default = 1)
% plt  = create PDF (=1) or PDF + randomized plot (= 2)            (default = 0)
% iter = # iterations for randomization test                       (default = 0)
%
% G      = Greenwood's statistic
% n      = sample size
% p      = parametric p-value
% p_rand = randomized p-value
%
% SEE ALSO: f_greenwood_cdf, f_greenwood_pdf, f_greenwood_rnd

% -----Notes:-----
% Greenwood's (1946) statistic (G) provides a relative measure of uniformity in
% the spacings among n ordered values that represent the occurrence of events in
% time or the placement of objects in space. Very small values of G indicate
% super-uniformity, while large values suggest departure from the random uniform
% spacing expected under the null hypothesis.
%
% In this function n is number of observations in X minus 2 and n+1 is the
% number of intervals in D that are squared & summed to produce G. Note that the
% sum of D = 1.
%
% For complete uniformity, the expected value of D_i = 1/(n+1), where n is the
% sample size and the corresponding value of G = 1/(n+1).

% -----References:-----
% D'Agostino, R. B. and M. A. Stephens. 1986.  Goodness of fit techniques.
%  Marcel Dekker, Inc. New York.
% David, H. A. and H. N. Nagaraja. 2003. Order Statistics. Wiley Series in
%  Probability and Statistics. John Wiley & Sons, Inc. 3rd Edition. (p. 133 in
%  Section 6.4: Random Division of an Interval).
% Greenwood, M. 1946. The statistical study of infectious diseases. J. R.
%  Statist. Soc. A 109(2): 85-110.
% Koziol, J. A. 1987. An alternative formulation of Neyman's smooth goodness of
%  fit test under composite alternatives. Metrika 34: 17-24.
% Lafaye de Micheaux, P. and V. A. Tran. 2014. PoweR: Reproducible research
%  tool to ease Monte-Carlo power simulation studies for goodness-of-fit tests
%  in R. Journal of Statistical Software.
% Stephens, M. A. 1981. Further Percentage Points for Greenwood's Statistic.
%  Journal of the Royal Statistical Society. Series A (General) 144(3): 364-366

% -----Implementation:-----
% Here's a couple of ways to implement the statistic:
%
% --METHOD 1 (follows David & Nagaraja, 2003 & Stephens, 1981):--
% Transform values to a 0-1 interval:
% U = unifcdf(X,min(X),max(X)); % after 'stat70.cpp' in PoweR package
%
% Get length of n+1 consecutive intervals between ordered values:
% D = diff([0;sort(U);1]);
%
% Sum of n+1 squared distances (Greenwood's Statistic):
% G = sum(D.^2);
%
% This is can be accomplished with 1 line of code:
% G(1) = sum(diff([0;sort(unifcdf(X,min(X),max(X)));1]).^2);
%
%
% --METHOD 2 modified for values that don't represent a fraction of the total:--
% (after D'Agostino & Stephens, 1986):
%
% Get length of consecutive intervals between ordered values:
% D = diff(sort(X));
% G = sum(D.^2) / sum(D).^2;
% 
% --METHOD 3: after (Askew, 1983):--
% X    = sort(X);                   % 1) ordered distances
% U    = (X./X(end)).^2;            % 2) transformed distances
% G(1) = sum([U(1)^2; diff(U).^2]); % 3) Greenwood's statistic
% 
% -> however, squaring the values in step 2 of Askew's method departs from the
%    actual PDF/CDF of Greenwood's statistic.

% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Apr-2014: now calculates parametric p-value; improved plotting
% Nov-2014: updated documentation

% -----Check input & set defaults:-----
if (nargin < 2), verb = 1; end % default don't send output to display
if (nargin < 3), plt  = 0; end % default don't create a plot
if (nargin < 4), iter = 0; end % default no randomization test

% Check input data:
[nr,nc] = size(X);
if (nc>1)
   error('X must be a column vector!');
end

if (nr<3), error('Raw data must have at least 3 rows!'); end

% Check plot options:
if (plt>1 && iter<0)
   error('ITER must be > 0 if PLT > 1');
end
% -------------------------------------

% Preallocate:
if (iter>0)
   G = nan(iter,1);
else
   G = NaN;
end

% Calculate Greenwood's statstic:
n    = nr-2;                  % sample size
D    = diff(sort(X));         % length of intervals between ordered observations
G(1) = sum(D.^2) / sum(D).^2; % observed statistic

if (nargout<3 && verb<1), return; end

% Test null hypothesis that data are homogeneous:
p = 1-f_greenwood_cdf(G(1),n);

% -----Optional Randomization Test:-----
if (iter>0)
   for i=2:(iter)
      R    = f_randWH(nr,0);        % uniform random numbers from 0-1
      D    = diff(sort(R));
      G(i) = sum(D.^2) / sum(D).^2; % randomized statistic
   end
   p_rand = sum(G>=G(1))/(iter);    % convert counts to probability
else
   p_rand = NaN;
end
% -----------------------------

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('              Greenwood''s Statistic:             \n');
   fprintf('--------------------------------------------------\n');
   fprintf('G stat  = %-3.4f (n = %d) \n',G(1),n);
   fprintf('p-value = %3.5f \n',p);
   if (iter>0)
      fprintf('p_rand  = %3.5f (with %d iterations)\n',p_rand,iter);
   end
   fprintf('--------------------------------------------------\n');
end

% -----Create parametric plot:-----
if (plt>0)
   % Plot Probability Distribution Function:
   [~,~,hdl(1)] = f_greenwood_plt(n,'PDF',1);
   
   % Change the title:
   title('');
   
   % Add observed statistic:
   hold on;
   hdl(2) = plot(G(1),0);
   set(hdl(2),'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
      'MarkerFaceColor','k','MarkerSize',8);
   legend(hdl,'PDF','Observed \itG','location','NorthEast')
end
% ------------------------

% -----Create a plot:-----
if (plt>1 && iter>0)
   
   colorVar = [1 1 1]*0.65; % histogram color
   
   figure; set(gcf,'color','w');
   hold on;
   box on;
   grid on;
   
   % Plot histogram:
   Y      = G(2:end);
   nbins  = max(min(length(Y)./10,100),50); % get # bins
   yi     = linspace(min(Y),max(Y),nbins)'; % range
   dy     = mean(diff(yi));                 % step size
   yfi    = histc(Y,yi-dy);                 % get frequencies
   yfi    = yfi./sum(yfi)./dy;              % scale frequencies
   hdl(1) = bar(yi,yfi,'FaceColor',colorVar,'EdgeColor','none');
   
   % Plot observed statistic:
   hdl(2) = plot(G(1),0);
   set(hdl(2),'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
      'MarkerFaceColor','k','MarkerSize',8);
   
   % Customize plot:
   xTxt = sprintf('Greenwood''s Statistic (n = %d)',n);
   xlabel(xTxt);
   xlabel(xTxt);
   ylabel('Empirical Probability Density');
   legend(hdl,'EDF','Observed \itG','location','NorthEast')
   minG = 1/(n+1);
   xlim([minG-(minG*0.1) 1]);
   % --------------------------
end

% Wrap results up into a structure:
result.G = G(1);           % oberved statistic
result.n = n;              % sample size
result.p = p;              % p-value
if (iter>0)
   result.p_rand = p_rand; % randomized p-value
end
