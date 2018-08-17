function [out,p] = f_outlyingness_test(X,method,iter,verb,plt)
% - permutation test of maximum outlyingness
% 
% USAGE: out = f_outlyingness_test(X,method,iter,verb,plt)
% 
% X      = input data (rows = observations, cols = variables)
% method = as: 'raw' or 'uni'
% iter   = # iterations for permutation test                       (default = 0)
% verb   = optionally send results to display                      (default = 0)
% plt    = optionally create a plot                                (default = 0)
% 
% out = maximum observed outlyingness measure
% p   = permutation-based p-value 

% -----Notes:-----


% -----References:-----
% Breiman, L., and A. Cutler. 2003. Manual on setting up, using, and
%  understanding Random Forests v4.0. Technical Report. 
%  ftp://ftp.stat.berkeley.edu/pub/users/breiman/Using_random_forests_v4.0.pdf

% -----Author:-----
% by David L. Jones, May-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), iter = 0; end % no randomization test by default
if (nargin < 4), verb = 0; end % default don't send output to display
if (nargin < 5), plt  = 0; end % default don't create a plot

% Force lower case:
method = lower(method);
% -------------------------------------

% Get # observatins:
n = size(X,1);

% Preallocate:
if (iter>0)
   out = nan(iter,1);
else
   out = NaN;
end

% Get observed statistic:
out(1) = max(f_outlyingness(X));

% Get smoothing parameter & plot title:
switch method
   case 'raw'
      h    = [];
      tTxt = ('Permutation method: raw data');
   case 'uni'
      h    = 1/sqrt(n); % Classell (2010)
      tTxt = sprintf('Perturbation method: %s',method);
   otherwise
      error('Unknown METHOD!')
end

% Perform significance test:
if (iter>0)
   switch method
      case 'raw' % Permute RAW data
         for i=2:iter
            out(i) = f_outlyingness(f_shuffle(X,3));
         end
      case 'uni' % Perturb univariate data
         for i=2:iter
            X_p    = X + randn(n,1)*h; % eq. 6.15 in Silverman (1986)
            out(i) = max(f_outlyingness(X_p));
         end
      otherwise
         error('Unknown METHOD!')
   end
   p = sum(out>=out(1))/(iter); % convert counts to probability
else
   p = NaN;
end

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('                 Outlier Detection:               \n');
   fprintf('--------------------------------------------------\n');
   fprintf('Outlyingness = %-3.4f (n = %d) \n',out(1),n);
   fprintf('p-value      = %3.5f \n',p);
   fprintf('# iterations = %d \n',iter);
   fprintf('method       = %s\n',method);
   if ~isempty(h)
      fprintf('h            = %3.5f\n',h);
   end
   fprintf('--------------------------------------------------\n');
end

% -----Create a plot:-----
if (plt>0 && iter>0)
   
   colorVar = [1 1 1]*0.65; % histogram color
   
   figure; set(gcf,'color','w');
   hold on;
   box on;
   grid on;
   
   % Plot histogram:
   Y      = out(2:end);
   nbins  = max(min(length(Y)./10,100),50); % get # bins
   yi     = linspace(min(Y),max(Y),nbins)'; % range
   dy     = mean(diff(yi));                 % step size
   yfi    = histc(Y,yi-dy);                 % get frequencies
   yfi    = yfi./sum(yfi)./dy;              % scale frequencies
   hdl(1) = bar(yi,yfi,'FaceColor',colorVar,'EdgeColor','none');
   
   % Plot observed statistic:
   hdl(2) = plot(out(1),0);
   set(hdl(2),'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
      'MarkerFaceColor','k','MarkerSize',8);
   
   % Customize plot:
   title(tTxt,'Interpreter','none');
   xTxt  = sprintf('Outlyingness Measure (n = %d)',n);
   xlabel(xTxt);
   xlabel(xTxt);
   ylabel('Empirical Probability Density');
   legend(hdl,'EDF','Observed Measure','location','NorthEast')
   % --------------------------
end

% Format for output:
out = out(1);