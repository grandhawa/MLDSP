function result = f_greenwood_rnd(n,iter,plt)
% - randomized distribution of Greenwood's statistic
%
% USAGE: result = f_greenwood_rnd(n,iter,plt);
%
% n    = sample size
% iter = # of randomization iterations                          (default = 1000)
% plt  = plot the frequency of G (=1) or H (=2)                    (default = 1)
%
% result = structure of results with the following fields:
%  .G  = Greenwood's Statistic
%  .H  = G scaled to 0-1 interval
%  .n  = sample size
%
% SEE ALSO: f_greenwood, f_greenwood_cdf, f_greenwood_pdf, f_randWH

% -----References:-----
% Moran, P. A. P. 1953. The random division of an interval--Part III. J. R.
%  Statisti. Soc. B 15: 77-80.

% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Apr-2014: N now correctly represents 'sample size'; improved setting of axis
%           limits

% -----Set defaults & check input:-----
if (nargin < 2), iter = 1000; end % default 1000 iterations
if (nargin < 3), plt  = 0;    end % default don't create a plot

% Check n:
if ~isscalar(n), error('N must be a scalar!'); end
if (n<1), error('N must be greater than 0'); end
if (n-floor(n)>0), error('N must be a whole number'); end
% ------------------------------------

nr = n+2; % get # rows

G = nan(iter,1); % preallocate
for i=1:(iter)
   R    = f_randWH(nr,0); % uniform random numbers from 0-1
   temp = f_greenwood(R,0);
   G(i) = temp;           % Greenwood's statistic
end

% Modified form of G that always ranges from 0-1:
H = ((n+1)/n)*(G - 1./(n+1)); % Moran (1953)

% -----Optionally create a plot-----
if (plt>0)
   
   switch plt
      case 1
         Y = G;
         xTxt = sprintf('Greenwood''s Statistic  (n = %d)',n);
      case 2
         Y = H;
         xTxt = sprintf('Index of Homogeneity  (n = %d)',n);
      otherwise
         error('Unknown PLT option!')
   end
   
   colorVar = [1 1 1]*0.65; % define histogram color
   
   figure; set(gcf,'color','w');
   hold on;
   box  on;
   grid on;
   
   % Plot histogram:
   nbins = max(min(length(Y)./10,100),50); % get # bins
   yi    = linspace(min(Y),max(Y),nbins)'; % range
   dy    = mean(diff(yi));                 % step size
   yfi   = histc(Y,yi-dy);                 % get frequencies
   yfi   = yfi./sum(yfi)./dy;              % scale frequencies
   hdl   = bar(yi,yfi,'FaceColor',colorVar,'EdgeColor','none');
   
   % Customize plot:
   title({'Randomized Distribution','of Uniform Spacing'});
   xlabel(xTxt);
   ylabel('Relative Frequency');
   
   % Set axis limits H:
   switch plt
      case 1
         minG       = 1/(n+1);
         axisVar    = axis;
         dx         = (1 - minG)* 0.05;
         axisVar(1) = minG - dx;
         
         if axisVar(2)>1, axisVar(2)= 1; end
         axis(axisVar);
      case 2
         axisVar = axis;
         axis([0 1 axisVar(3:4)])
   end
   
end
% ----------------------------------

% Wrap results up into a structure:
result.G = G;
result.H = H;
result.n = n;
