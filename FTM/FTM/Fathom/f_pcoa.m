function result = f_pcoa(dis,diag,scale,neg)
% - Principal Coordinates Analysis (PCoA)
%
% USAGE: result = f_pcoa(dis,diag,scale,neg);
%
% dis    = square symmetric distance matrix.
% diag   = make diagnostic plots                                   (default = 0)
% scale  = scale eigenvectors (= scores) by their eigenvalue       (default = 1)
% neg    = discard (= 0), keep (= 1), or correct (= 2)
%          negative eigenvalues                                     (default = 0)
% 
% result = structure of results with the following fields:
%  .scores = scaled (default) or normalized eigenvectors
%  .evals  = eigenvalues
%  .expl   = percent and cumulative variance explained
%
% SEE ALSO: f_pcoaPlot, f_nmds, f_pca

% -----References:-----
% Anderson, M. J. 2002. PCOORD: a FORTRAN computer program for principal
%   coordinate analysis. Dept. of Statistics University of Auckland.
%   Available from: http://www.stat.auckland.ac.nz/PEOPLE/marti/
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.
% Legendre, P., and M. J. Anderson. 1999. Distance-based redundancy analysis:
%   testing multispecies responses in multifactorial ecological experiments.
%   Ecol. Monogr. 69: 1?24.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.

% -----Notes:-----
% Principal coordinates associated with negative eigenvalues are now optionally
% corrected using Method 1 of Legendre & Legendre (1998). Note that negative
% eigenvalues arise when attempting to use PCoA to obtain a Euclidean embedding
% of a semi-metric dissimilarity coefficients (e.g., Bray-Curtis) and represent
% negative sums-of-squares (SS) needed to offset the inflated SS of the positive
% eigenvalues. Thus, the SS Total will be inflated if (1) you retain only the
% positive eigenvalues or (2) you correct for negative eigenvalues by adding a
% constant as in the method of Legendre & Legendre (1998). See McArdle &
% Anderson (2001) for more details.

% -----Author:-----
% by David L. Jones, Apr-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: extended plots are now optional
% Mar-2008: num2Str changed to num2str
% Apr-2008: calculate scores after (vs. before) optionally remove neg
%           eigenvalues; updated doc; set tolerance for 0 eigenvalues
% Mar-2010: added support for optionally correcting for negative eigenvalues
% May-2012: results returned as a structure; plotting routines moved to a
%           separate function (f_pcoaPlot)
% Dec-2012: updated documentation

% -----Check input & set defaults:-----
if (nargin < 2), diag   = 0; end % default don't create diagnostic plots
if (nargin < 3), scale  = 1; end % default scale eigenvectors
if (nargin < 4), neg    = 0; end % default discard negative eigenvalues


% Check if symmetric distance matrix:
if (f_issymdis(dis) == 0)
   error('Input DIS must be square symmetric matrix');
end

% Set tolerance for eigenvalues > 0:
tol = sqrt(eps);
% -------------------------------------

n = size(dis,1);

% Gower's centered matrix:
uno = ones(n,1);
I   = eye(n,n);
A   = -0.5*(dis.^2);
G   = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno')); % Anderson, 2002

% -----Eigenanalysis:-----
[evects,evals] = f_eig(G);

% Discard eigenvalues = 0:
idx    = find (abs(evals)>tol);
evects = evects(:,idx);
evals  = evals(idx);

% Only need n-1 axes for n objects:
if (size(evals,1) > n-1)
   evects = evects(:,1:(n-1));
   evals  = evals(1:(n-1));
end

% Variance explained (same as Anderson's PCOORD):
varExplained  = (evals/sum(evals))*100;       % Percent Variance Explained
cvarExplained = cumsum(varExplained);         % Cumulative Percent Explained
expl          = [varExplained cvarExplained]; % combine into 1 matrix

% Handle negative eigenvalues:
nNeg = sum(find(evals<0)); % get # of neg eigenvalues
if (nNeg>0) && (neg==2)    % correct
   c                   = abs(min(evals));
   disCOR              = f_rewrap(sqrt( f_unwrap(dis).^2 + 2*c )); % eq. 9.25 L&L, 1998
   temp   = f_pcoa(disCOR,0,0,1);% redo PCoA, don't scale yet
   % Extract from structure:
   evects = temp.scores;
   evals  = temp.evals;
   expl   = temp.expl;
   clear temp;
elseif (nNeg>0) && (neg<1) % discard
   idx    = find(evals>0);
   evects = evects(:,idx);
   evals  = evals(idx);
   expl   = expl(idx,:);   
end

% Scale eigenvectors (= scores):
if (scale>0)
   % Scale eigenvectors by their eigenvalue to mean = 0,
   % stdv = eigenvalue/(n-1), and length = eigenvalue/(n-1):
   scores = evects .* repmat(abs((evals.^0.5)'),size(evects,1),1);
else
   % Normalized eigenvectors (lengths = 1):
   scores = evects;
end

% -----Create optional diagnostic plots:-----
if (diag>0) 
   % Shepard Plot of 2-d configuration:
   figure;
   set(gcf,'color','w'); % set bg color to white
   % Euclidean distance matrix 1st 2 eigenvectors as a vector:
   edist     = f_dis(real(scores(:,1:2)),'euc'); % discard imaginary portions
   Rs        = f_corr(f_unwrap(dis),f_unwrap(edist))^2;
   titleVar  = ['2-d Shepard Diagram (R^2 = ' num2str(Rs) ')'];
   xlabelVar = ['Original Dissimilarites (' num2str(n) ' objects)'];
   ylabelVar = 'Fitted Distances (1st two eigenvectors)';
   plot(dis, edist, 'bo');
   title(titleVar);
   xlabel(xlabelVar);
   ylabel(ylabelVar);
   grid on;
   
   % Scree Plot of residual variance:
   figure;
   set(gcf,'color','w');               % set bg color to white
   stem(100-expl(:,2),'-r.');hold on;  % plot residual variance
   plot(100-expl(:,2),'-b'); hold off; % ditto
   axis([0.25 (size(expl,1)+ 0.75) 0 1.08*max(100-expl(:,2))]);
   title('Scree Plot: Variance Explained vs. Dimensionality');
   xlabel('# Dimensions');
   ylabel('Residual Variance (%)');
   set(gca,'Xtick',1:(size(expl,1)));
end
% -------------------------------------------

% Wrap results up into a structure:
result.scores = scores; % scaled (default) or normalized eigenvectors
result.evals  = evals;  % evals  = eigenvalues
result.expl   = expl;   % percent and cumulative variance explained
