function [G,den,hdl] = f_greenwood_plt(n,type,plt)
% - plot the PDF or CDF for Greenwood's statistic
% 
% USAGE: [G,den,hdl] = f_greenwood_plt(n,type,plt);
% 
% n    = sample size for Greenwood's statistic
% type = type of plot to create as: 'PDF' or 'CDF'
% plt  = create a plot                                             (default = 1)
% 
% G   = Greenwood's statistic
% den = (cumulative) probability densities
% hdl = graphics handle
% 
% SEE ALSO: f_greenwood, f_greenwood_cdf, f_greenwood_pdf

% -----Notes:-----
% The PLT option allows you to calculate the densities of G for a given N
% without plotting then.

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), plt = 1; end % default create a plot

% Check N:
if ~isscalar(n), error('N must be a scalar!'); end
if (n<1), error('N must be 1 or greater!'); end
if (n-floor(n)>0), error('N must be a whole number'); end
% -------------------------------------

type     = upper(type); % force upper case
nPts     = 1000;        % set # points
widthVar = 1;           % define line width

% Get range of G statistic:
G = linspace(1/(n+1),1,nPts)';

% Calculate probability densities:
switch type
   case 'PDF'
      den = f_greenwood_pdf(G,n);
      tTxt = 'PDF';
      yTxt = 'Probability Density';
   case 'CDF'
      den = f_greenwood_cdf(G,n);
      tTxt = 'CDF';
      yTxt = 'Cumulative Probability Density';
   otherwise
      error('Unknown plot type!');
end

% Create plot:
if (plt>0)
   figure; set(gcf,'color','w');
   hdl = plot(G,den,'k-','LineWidth',widthVar);
   title(tTxt);
   xTxt = sprintf('  (n=%d)',n);
   xlabel(['Greenwood''s Statistic' xTxt]);
   ylabel(yTxt);
   grid on; box on;
   
   % Set axis limits:
   minG = 1/(n+1);
   xlim([minG-(minG*0.1) 1]);
end
