function result = f_variogram(yDis,X,n,iter,plt)
% - multivariate empirical variogram
%
% USAGE: result = f_variogram(yDis,X,n,iter,plt);
%
% yDis = square symmetric distance matrix derived from response variables
% X    = matrix of site coordinates
% n    = # of distance classes, use n=0 for Sturges Rule           (default = 0)
% iter = # iterations for permutation test                         (default = 0)
% plt  = create plot                                               (default = 0)
% 
% result = structure of results with the following fields:
%  .dxy  = sequential bounds of distance classes
%  .d    = center of each distance class
%  .C    = coefficient of complementarity for each distance class
%  .p    = boolean indicating which values of 'C' are significant (=1)
%  .cv   = critical values for permutation-based confidence interval
%  .tot  = total inertia of the response variable
%  .nc   = number of connections within each distance class
%  .nw   = number of members within each distance class
% 
% SEE ALSO: f_correlogram, f_selectW, f_moran

% -----Notes:-----
% Variograms decompose the spatial variability of a variable among distance
% classes, providing a non-standardized version of Geary's coefficient of
% spatial autocorrelation (Legendre & Legendre, 1998). This function creates a
% variogram of complementarity (C), which provides a coefficient of the
% dissimilarity in species composition among sites that are members of the same
% distance class (Wagner, 2003).
% 
% Total inertia of the response variable (tot) provides a baseline that
% represents the average squared distance among pairs of sites. Values of C <
% total inertia indicates that the average squared distance between sites within
% a distance class is LESS the average distance among all sites (= positive
% autocorrelation). Values of C > total inertia indicates that the average
% squared distance between sites within a distance class is MORE the average
% distance among all sites (= negative autocorrelation). Values of C = total
% inertia indicates NO correlation.
% 
% Note that due to bias, variograms are not generally interpreted beyond half
% the maximum distance between sites, which is indicated by a dotted line in the
% variogram plot (Wagner, 2003; Dray et al., 2006).

% -----Reference:-----
% Borcard, D., F. Gillet, and P. Legendre. 2011. Numerical Ecology with R.
%  Springer, NY.
% Dray, S. 2006. Moran's eigenvectors of spatial weighting matrices in R.
%   Documentation included in the 'spacemakeR for   R' documentation. Available
%   from: http://biomserv.univ-lyon1.fr/~dray/software.php
% Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial-modelling: a
%  comprehensive framework for principal coordinate analysis of neighbor matrices
%  (PCNM). Ecological Modelling 196: 483-493.  
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% Wagner, H. H. 2003. Spatial covariance in plant communities: integrating
%  ordination, geostatistics, and variance testing. Ecology 84(4): 1045-1057.

% -----Author:-----
% by David L. Jones, Mar-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 2011-Mar: replaced f_euclid with f_dis, support for Sturges rule, input
%           now requires response variable to be a distance matrix

% -----Set defaults & check input:-----
if (nargin < 3), n    = 0; end % default use Sturges Rule
if (nargin < 4), iter = 0; end % default no permutation test
if (nargin < 5), plt  = 0; end % default no plot

% Check if yDis & X are compatible:
if size(yDis,1) ~= size(X,1), error('yDis & X need same # of rows'); end

% Check imput if this isn't a permutation run:
if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix');
end
% -------------------------------------

xDis = f_dis(X,'euc'); % Euclidean distance of spatial variable
dMax = max(xDis(:));   % max distance
dMin = 0;              % default min distance
nr   = size(yDis,1);   % # rows in response variable

% Sturges Rule (Borcard et al. 2011: p.237):
if (n==0)
   n = ceil(1 + (3.3219 * log10(nr*(nr-1)/2)));
end

% Total inertia of the response variable:
dis     = yDis.^2;
SSt_tot = sum(dis(:))/2;
tot     = SSt_tot / (nr*(nr-1)/2) / 2; % Wagner, 2003: paragraph after eq. 7

% Construct distance classes:
dxy = linspace(dMin,dMax,n+1)';

% Put NaN's along the diagonal of the symmetric distance matrix, so sites 
% paired with themselves are not considered members of the distance class = 0
% when calculating 'nw':
xDis(eye(size(xDis))==1) = NaN;

% Complementarity for each distance class:
bin{n} = zeros(size(xDis)); % preallocate
SSt    = zeros(n,1);
nc     = zeros(n,1);
nw     = zeros(n,1);
d      = zeros(n,1);
C      = zeros(n,1);
cv     = zeros(n,2);
p      = zeros(n,1);
% 
for i = 1:n
   % Get centers of each distance class:
   d(i) = mean([dxy(i) dxy(i+1)]);

   % Binary symmetric matrix specifying membership in each distance class:
   bin{i} = ( xDis>=dxy(i) & xDis<=dxy(i+1) );

   % Sum-of-Squares total within each distance class (~ eq. 13.10; L&L, 1998):
   dis    = (yDis.^2 .* bin{i});
   SSt(i) = sum(dis(:))/2;

   % Number of connections within each distance class:
   nc(i) = sum(sum(bin{i},2) > 0);

   % Number of members within each distance class:
   nw(i) = (sum(bin{i}(:)))/2;

   % Complementarity:
   C(i) = SSt(i)/nw(i)/2;
end

% Get half the maximum distance between sites:
siteDis = xDis(:);
idx     = find(siteDis == max(siteDis));
hMaxDis = siteDis(idx(1))/2;


% -----Permutation Test:-----
if (iter>0)
   fprintf('Permuting the data %d times...\n',iter-1);
   Cperm = zeros(n,iter-1);      % preallocate
      
   for i = 1:(iter-1)                    % observed value is considered a permutation
      yDisPerm      = f_shuffle(yDis,2); % permute observations
      resultPerm = f_variogram(yDisPerm,X,n,0,0);
      Cperm(:,i) = resultPerm.C;         % permuted C
   end

   % Specify the significance level:
   alpha = 0.05;
   
   for k = 1:n
      % Critical values for confidence interval (Bonferroni-corrected):
      cv(k,1) = prctile(Cperm(k,:),((alpha/2)/n)*100);
      cv(k,2) = prctile(Cperm(k,:),(1-(alpha/2)/n)*100);
      
      % Determine if C is significant:
      p(k) = logical(C(k)<=cv(k,1)  || C(k)>=cv(k,2));
   end
else
   p  = repmat(NaN,n,1);
   cv = repmat(NaN,n,2);
end
% ---------------------------


% -----Wrap results up into a structure:-----
result.dxy = dxy;
result.d   = d;
result.C   = C;
result.p   = p;
result.cv  = cv;
result.tot = tot;
result.nc  = nc;
result.nw  = nw;


% -----Create plot:-----
if (plt>0)
   figure;
   set(gcf,'color','w'); % set background color to white
   hold on;
   
   if (iter>0) % Permutation test:
      plot(d,C,'b-');                             % Distance vs. Complementarity
      h(1) = plot(d(p==1),C(p==1),'b.','MarkerSize',18); % significant C's
      h(2) = plot(d(p==0),C((p==0)),'bo');               % non-significant C's
      h(3) = plot(d,cv(:,1),'k-');                       % confidence interval
             plot(d,cv(:,2),'k-');

      axisVar = axis;
      h(4) = plot([hMaxDis;hMaxDis],[axisVar(3);axisVar(4)],'k:'); % half max dist
      h(5) = plot([axisVar(1);axisVar(2)],[tot;tot],'k--');        % total inertia

      legend(h,'Significant','Non-significant','Confidence Interval',...
         'Half Max Distance','Total Inertia');

   else % No permutation test:
      h(1) = plot(d,C,'bo-'); % Distance vs. Complementarity
      axisVar = axis;
      h(2) = plot([hMaxDis;hMaxDis],[axisVar(3);axisVar(4)],'k:'); % half max dist
      h(3) = plot([axisVar(1);axisVar(2)],[tot;tot],'k--');      % total inertia

      legend(h, 'Variogram', 'Half Max Distance','Total Inertia');
   end

   xlabel('\bfDistance')
   ylabel('\bfComplementarity');
   box on;
end

