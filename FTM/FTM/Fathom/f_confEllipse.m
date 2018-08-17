function [h,exy] = f_confEllipse(xy,conf,iter,plt,cSpec)
% - parametric or bootstrapped confidence ellipse for bivariate data
%
% USAGE: f_confEllipse(xy,conf,iter,plt,colorVar)
%
% xy    = 2 column matrix of input data
% conf  = confidence interval enclosed by ellipse (default = 0.95)
% iter  = # iterations for bootstrap resampling   (default = 0)
% plt   = add ellipse to current plot             (default = 0)
% cSpec = line spec defining ellipse color        (default = 'r');
% 
% h   = handle to plotted ellipse
% exy = coordinates of ellipse
%
% SEE ALSO: f_convHull

% -----Notes:-----
% if ITER = 0 a parametric confidence ellipse is generated, otherwise a
% bootstrapped version is computed. For bootstrapped ellipses it is
% recommened that ITER be at least 1000.

% -----References:-----
% modified after Douglas M. Schwarz's CONFELLIPSE2 
% posted to news://comp.soft-sys.matlab, 1998

% -----Author:-----
% by David L. Jones, Dec-2003
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% May-2008: replaced | with ||
% Aug-2009: preallocate Fstat, updated call to f_mregress which now returns a
%           structure
% Oct-2010: updated documentation

% -----Set defaults and check input:-----
if (nargin < 2), conf  = 0.95; end; % default 95% confidence ellipse
if (nargin < 3), iter  = 0;    end; % default 0 bootstrap samples
if (nargin < 4), plt   = 0;    end; % default don't plt ellipse
if (nargin < 5), cSpec = 0;    end; % default color 'r'

if (size(xy,2) ~= 2)
   error('XY must be a 2 column matrix')
end

if (conf <= 0) || (conf >= 1)
   error('CI must be >0 and <1)');
end
% ---------------------------------------

n   = size(xy,1);
mxy = mean(xy);

numPts = 181; % The number of points in the ellipse.
th     = linspace(0,2*pi,numPts)';
p      = 2;   % 2-Dimensional data set


% -----Get critical values for F-stat:-----
if (iter<1) % Parametric:
   
   cv = finv(conf,p,n-p);  
   
else % Bootstrap:
   
   % Extract XY components:
   x = xy(:,1);
   y = xy(:,2);
   
   Fstat = zeros(iter,1); % preallocate result array
      
   % Get observed F-statistic
   model    = f_mregress(x,y,0,0,0);
   Fstat(1) = model.F_stat;
   
   % Get bootstrapped F-statistic:
   for i=2:iter
      model    = f_mregress(x,f_boot(y),0,0,0);
      Fstat(i) = model.F_stat;
   end
   
   % Bootstrapped critical value:
   cv = prctile(Fstat,[1-(1 - conf)]*100);    
end
% ----------------------------

% Principal Axes:
k              = cv*p*(n-1)/(n-p);
[pc,score,lat] = princomp(xy);
ab             = diag(sqrt(k*lat));
exy            = [cos(th),sin(th)]*ab*pc' + repmat(mxy,numPts,1);

% Add ellipse to current plot:
if (plt>0)
   hold on;
   h = line(exy(:,1),exy(:,2),'Clipping','off');
   set(h,'Color',cSpec);
else
   h = NaN;
end




