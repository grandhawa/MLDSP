function result = f_correlogram(yDis,X,n,iter,mc,prog,trim,verb,plt)
% - multivariate Mantel correlogram
%
% USAGE: result = f_correlogram(yDis,X,n,iter,mc,prog,trim,verb,plt);
%
% yDis = square symmetric distance matrix derived from response variables
% X    = matrix of site coordinates
% n    = # of distance classes; use n=0 for Sturges Rule           (default = 0)
% iter = # iterations for permutation test                         (default = 0)
% mc   = correction used for multiple testing: 'holm' (= Holmes),
%        'bon' (= Bonferroni), 'ds' (= Dunn-Sidak),  or 'none'  (default = none)
% prog = progressive correction for multiple testing               (default = 1)
% trim = limit distance classes beyond half maximum distance       (default = 1)
% verb = send output to display                                    (default = 1)
% plt  = create plot                                               (default = 0)
%
% result = structure of results with the following fields:
%  .dxy  = sequential bounds of distance classes
%  .d    = center of each distance class
%  .C    = coefficient of complementarity for each distance class
%  .p    = boolean indicating which values of 'C' are significant (= 1)
%  .cv   = critical values for permutation-based confidence interval
%  .tot  = total inertia of the response variable
%  .nc   = number of connections within each distance class
%  .nw   = number of members within each distance class
%
% SEE ALSO: f_adjustP, f_variogram, f_mantel

% -----Notes:-----
% In this function, positive Mantel statistics (+r) correspond to
% positive spatial autocorrelation, while negative Mantel statistics (-r)
% correspond to negative spatial autocorrelation.
% 
% Note that due to bias, correlograms are generally not interpreted beyond half
% the maximum distance between sites, which is indicated by a dotted line in the
% variogram plot (Wagner, 2003; Dray et al., 2006; Borcard et al., 2011).
% 
% TRIM=1: limits distance classes beyond half the maximum distance to only
% those that include all points; TRIM=0 includes all distance classes.
%
% When PROG=1, a progressive correction for multiple comparison tests is
% used following the method outlined in Legendre & Legendre (1998: p.721).
% Note they recommend using a 'progressive' correction when it is known
% there is significant spatial autocorrelation present in the first
% (smallest) distance class (i.e., no correction is applied to the first
% p-value) and one needs to determine which distance class the
% autocorrelation extends.
%
% This function has been tested against the mantel.correlog.R function from
% the vegan package for R and provides similar results.

% -----Reference:-----
% Borcard, D., F. Gillet, and P. Legendre. 2011. Numerical Ecology with R.
%   Use R! series. Springer, NY.
% Borcard, D., and P. Legendre. 2012. Is the Mantel correlogram powerful
%   enough to be useful in ecological analysis? A simulation study. Ecology
%   93: 1473?1481.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.(736-738)
%
% mofified after f_variogram and mantel.correlog.R (the latter by P.
% Legendre in the vegan package for R)

% -----Author:-----
% by David L. Jones, Apr-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), n    = 0; end % default use Sturges Rule
if (nargin < 4), iter = 0; end % default no permutation test
if (nargin < 5), mc   = 0; end % default use Holmes correction method
if (nargin < 6), prog = 1; end % default progressive correction
if (nargin < 7), trim = 1; end % default distance classes to half maximum
if (nargin < 8), verb = 1; end % default send output to display
if (nargin < 9), plt  = 0; end % default create plot

% Check if yDis & X are compatible:
if size(yDis,1) ~= size(X,1), error('yDis & X need same # of rows'); end

% Check input:
if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix');
end
% -------------------------------------

xDis = f_dis(X,'euc'); % Euclidean distance of spatial variable
vDis = f_unwrap(xDis); % unwrap distance matrix as a vector
dMax = max(vDis);      % max distance
dMin = min(vDis);      % min distance
nr   = size(X,1);      % # obserations

% Sturges Rule (Borcard et al. 2011: p.237):
if (n==0)
   n = ceil(1 + (3.3219 * log10(nr*(nr-1)/2) ));
end

% Construct distance classes:
dxy = linspace(dMin,dMax,n+1)';

% Put NaN's along the diagonal of the symmetric distance matrix, so sites
% paired with themselves are not considered members of the distance class = 0
% when calculating 'nw':
xDis(eye(size(xDis))==1) = NaN;

% Get half the maximum distance between sites:
siteDis = xDis(:);
idx     = find(siteDis == max(siteDis));
hMaxDis = siteDis(idx(1))/2;

% Compute Mantel statistics for each distance class:
d      = zeros(n,1);
bin{n} = zeros(size(xDis)); % preallocate
nc     = zeros(n,1);
nw     = zeros(n,1);
r      = NaN(n,1);
p      = NaN(n,1);

for i = 1:n
   % Get centers of each distance class:
   d(i) = mean([dxy(i) dxy(i+1)]);
        
   % Binary symmetric matrix specifying membership in each distance class
   bin{i} = ( xDis>=dxy(i) & xDis<=dxy(i+1) );
   
   % Optionally trim distance classes > half max distance
   if (trim>0)
      if (d(i)> hMaxDis) && (sum(sum(bin{i}))==0)
         break
      end
   end
      
   % Number of connections within each distance class:
   nc(i) = sum(sum(bin{i},2) > 0);
   
   % Number of members within each distance class:
   nw(i) = (sum(bin{i}(:)))/2;
   
   % Mantel Test:
   if (iter>0)
      fprintf('Permuting the data %d times for %d of %d distance classes...\n',iter,i,n);
   end
   [r(i),p(i)] = f_mantel(yDis,bin{i},0,iter);
   
   % Reverse the sign of the Mantel Statistic (in this case 0 in the model
   % matrix indicates sites that belong to the SAME distance class; see Borcard
   % et al., 2011: p. 234 or Legendre & Legendre, 1998: p. 738 or Borcard &
   % Legendre, 2012: p1474):
   r(i) = r(i) * -1;
end

% Adjust p-values for multiple comparisons:
if (iter>0)
   if (prog>0) % Progressive correction:
      pA = NaN(n,1);
      pA(1) = p(1); % don't correct first value (L&L, 1998 p. 721)
      for j=2:n
         temp  = f_adjustP(p(1:j),mc);
         pA(j) = temp(j);
      end
   else % Simultaneous correction:
      pA = f_adjustP(p,mc);
   end
else
   pA = NaN(n,1);
end

% -----Create plot:-----
if (plt>0)
   
   alpha = 0.05; % specify the significance level
   
   figure;
   set(gcf,'color','w'); % set background color to white
   hold on;
   
   if (iter>0) % Permutation test:
      idx1 = pA<=alpha; % index to significant Mantel stats
      idx0 = pA>alpha;  % index to non-signifiant Mantel stats
      
      plot(d,r,'b-');                                       % Distance vs. Mantel Stat
      if (sum(idx1)>0)
         h(1) = plot(d(idx1),r(idx1),'b.','MarkerSize',18); % significant r's
      end
      h(2)    = plot(d(idx0),r((idx0)),'bo');               % non-significant r's
      axisVar = axis;
      h(3)    = plot([hMaxDis;hMaxDis],[axisVar(3);axisVar(4)],'k:'); % half max dist
      
      if (sum(idx1>0))
         legend(h,'Significant','Non-significant','Half Max Distance');
      else
         legend(h(2:3),'Non-significant','Half Max Distance');
      end
      
   else % No permutation test:
      h(1)    = plot(d,r,'bo-'); % Distance vs. Mantel Stat
      axisVar = axis;
      h(2)    = plot([hMaxDis;hMaxDis],[axisVar(3);axisVar(4)],'k:'); % half max dist
      
      legend(h, 'Correlogram', 'Half Max Distance');
   end
   
   % Customize plot:
   box on;
   xlabel('Distance Class')
   ylabel('Mantel Statistic');
   f_origin('h','-');
   title('\bfMultivariate Mantel Correlogram');
end

% -----Send output to display:-----
if (verb>0)
   fprintf('\n============================\n');
   fprintf(  '    Mantel Correlogram:\n');
   fprintf(  '============================\n');
   fprintf(  'Class:  r:     p:     pA:   \n');
   for i=1:n
      fprintf('%6.0f %+6.4f %6.4f %6.4f \n',i,r(i),p(i),pA(i))
   end
end
fprintf('\n');
fprintf('# iterations       = %d\n',iter);
fprintf('Correction Method = %s\n',mc)
fprintf('Progressive       = %d\n',prog);
fprintf(  '----------------------------\n');
fprintf('Class = distance class (half = %d) \n',hMaxDis);
fprintf('r     = Mantel statistic \n');
fprintf('p     = permutation-based p-value \n');
fprintf('pA    = adjusted p-value \n');
% ---------------------------------


% -----Wrap results up into a structure:-----
result.n    = n;   % # distance classes
result.dxy  = dxy; % distance classes
result.d    = d;   % center of each distance class
result.nc   = nc;  % # connections within each distance class
result.nw   = nw;  % # members of each distance class
result.r    = r;   % Mantel statistic
result.p    = p;   % corresponding permutation-based signifiance
result.pA   = pA;  % p-values adjusted for multiple comparison tests
