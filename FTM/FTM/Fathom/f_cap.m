function result = f_cap(Y,dis,X,W,sm,iter,verb,m,loo,b,opt)
% - Canonical Analysis of Principal Coordinates (CAP) using ANY distance matrix
%
% USAGE: result = f_cap(Y,'dis',X,W,sm,iter,verb,m,loo);
%
% Y   = matrix of response variables (rows = obs, cols = variables)
% dis = dissimilarity measure to apply to Y
%       (e.g., dis = 'bc'; see help for f_dis)
%
% X = (1) vector of integers specifying group membership for objects in yDis,
%     (2) ANOVA design matrix specified by dummy coding, or
%     (3) matrix of explanatory variables (rows = observations, cols = variables)
%
% W    = unknown observations to classify            (default = empty)
% sm   = use spatial median instead of centroid      (default = 0)
% iter = # iterations for permutation test           (default = 0)
% verb = optionally send results to display          (default = 1)
% m    = # PCoA axes to retain                       (default = 0)
% loo  = perform leave-one-out cross-validation      (default = 0)
%
% results = structure with the following fields:
%  .grp          = grouping vector
%  .Qm           = retained Principal Coordinates
%  .Qexpl        = variation of yDis explained by each PCoA
%  .siteScores   = coordinates in canonical space (= Qstar)
%  .fitScores    = fitted site scores
%  .resScores    = residal site scores
%  .centroids    = group centroids (or spatial medians) defined in X
%  .yCor         = correlation of canonical axes with original Y variables
%  .xCor         = correlation of canonical axes with X variables
%  .U            = canonical eigenvectors
%  .evals        = canonical eigenvalues (= canonical correlations)
%                  1st value is greatest root statistic)
%  .Ures         = residual eigenvectos
%  .evalsRes     = residual eigenvalues
%  .trc          = trace statistic
%  .grs          = greatest root statistic (1st squared canonical eigenvalue)
%  .p            = structure of randomized probabilities: p.trace (trace) and
%                  p.grs (greatest root)
%  .canVar       = fraction of variation explained by canonical axes
%  .resVar       = residual variance
%  .wQm          = predicted PCoA's of W
%  .RSS          = residual sum-of-squares for W (needed for loo_RSS)
%  .siteScores_W = site scores of W
%  .grpW         = predicted group membership of W
%  .PP           = posterior probabilities of group prediction
%  .marg         = measure of the confidence of classification of unknowns,
%                  values range from 0 (= lack) to 1 (= complete confidence)
%  .mDisW        = min distance-to-centroid doesn't exceed max distance within group (= 1)
%  .loo_err      = leave-one-out classification error rates
%  .loo_RSS      = leave-one-out residual sum-of-squares
%  .ctr          = central tendency ('centroid', 'spatial median', or 'none')
%  .method       = 'CDA' or 'CCorA'
%
% SEE ALSO: f_capOptimal, f_capPlot, f_chanceClass, f_dis, f_designMatrix

% -----Notes:-----
% This program performs nonparametric Canonical Discriminant Analysis (CDA) on
% ANY symmetric distance matrix when the input for X is (1) a vector specifying
% group membership or (2) an ANOVA design matrix. It performs nonparametric
% Canonical Correlation Analysis (CCorA) when X is (3) a matrix of explanatory
% variables.
%
% CDA finds canonical axes that maximizes the separation of groups specified
% in X, while CCorA finds axes that maximize the correlations with the
% explanatory variables in matrix X.
%
% You may use f_designMatrix to create an ANOVA design matrix for input as
% X; the matrix should be full rank (not singular) and DO NOT include an
% intercept term (a column of 1's). This is done for you if X is merely a
% vector of integers specifying group membership.
%
% The program asks the user to specify how many axes of Q to retain for the
% analysis (m). Examine the EIGENVALUES and '% VARIATION EXPLAINED' output in
% the command window and try to include as much information in Q as possible
% with as few axes as possible. 'Set LOO = 1' to and try different values
% of m to determine which results in the classification error (LOO_ERR) for
% CDA and/or the lowest residual sums-of-squares (LOO_RSS) for CCorA .
%
% The randomized probabilities given by P allow you to test the hypothesis of (1)
% no signficant differences among groups when X specifies group membership (CDA)
% or (2) no significant relationship with X when it is a matrix of explanatory
% variables (CCorA).
%
% Note: use siteScores to assess correlation of repsonse variables, and
% fitScores to assess correlation of explanatory variables.
%
% Note: 'b' is only passed internally by the function when performing a
% Leave-One-Out estimate of Residual Sums-of-Squares.
%
% Note: m < 0 when f_cap is called by f_capOptimal to determine the range of
% values of m to try.
% 
% Note: opt = 1 when f_cap is called by f_capOptimal or during the LOO

% -----References:-----
% Anderson, M. J. 2002. CAP: a FORTRAN program for canonical analysis of
%   principal coordinates. Dept. of Statistics University of Auckland. Available
%   from: http://www.stat.auckland.ac.nz/PEOPLE/marti/
% Anderson, M. J. & J. Robinson. 2003. Generalized discriminant analysis
%   based on distances. Aust. N. Z. J. Stat. 45(3): 301-318.
% Anderson, M. J. & T. J. Willis. 2003. Canonical analysis of principal
%   coordinates: a useful method of constrained ordination for ecology. Ecology
%   84(2): 511-525.
% Breiman, L., and A. Cutler. 2003. Manual on setting up, using, and
%   understanding Random Forests v4.0. Technical Report. Available from:
%   http://www.stat.berkeley.edu/~breiman/Using_random_forests_v4.0.pdf
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Jul-2002: Renamed f_cap (was f_cva) for consistency with Anderson's CAP;
%           Improved display of results
% Apr-2003: Overhauled code: f_eig-based, removed dependency on ortha.m,
%           follows Appendix of Anderson & Willis (2003) more closely.
% Mar-2008: Changed | to ||, & to &&, results returned as a structure
% Aug-2009: Changed num2Str to num2str, p returns NaN when iter=0
% Feb-2011: overhaul to closely follow the equations in Anderson &
%           Robinson (2003) and Anderson & Willis (2003); LOO CV is now
%           done correctly; RSS is now calculated; plotting moved to f_capPlot
% Mar-2011: removed SSt, SSr, SSe, and R2 as these apply to asymmetric
%           analyses like RDA
% Apr-2011: no longer returns Qexpl when m<0; make sure m is not > noAxes
% Feb-2013: set b = empty by default
% Mar-2013: use f_recode to convert GRP to consecutive integers; this is
%           critical to get high classification success rates!!
% Aug-2013: added an internal 'opt' flag so yCor and xCor are not calculated
%           when f_cap is called from f_capOptimal or during an LOO run.

% -----Set defaults & check input:-----
if (nargin <  4), W    = []; end % no observations to classify
if (nargin <  5), sm   = 0;  end % default use centroid vs. spatial median
if (nargin <  6), iter = 0;  end % no permutation test by default
if (nargin <  7), verb = 1;  end % send output to display by default
if (nargin <  8), m    = 0;  end % get m interactively
if (nargin <  9), loo  = 0;  end % no leave-one-out cross-validation by default
if (nargin < 10), b    = []; end % set b to empty by default
if (nargin < 11), opt  = 0;  end % set opt to zero by default

[n,ncY] = size(Y); % # obs, # variables

if n ~= size(X,1), error('Y & X need same # of rows'); end

% When X is a vector specifying group membership:
if (size(X,1)==1) || (size(X,2)==1) % row or col vector
   cda     = 1;                     % CDA:
   grp     = X(:);                  % make sure it's a column vector
   grp_old = grp;                   % keep a copy of original
   grp     = f_recode(grp);         % recode values as consecutive integers
   X       = f_designMatrix(grp,1); % create ANOVA design matrix w/o intercept term
   uGrp    = f_unique(grp);         % unique groups, unsorted
   nGrp    = length(uGrp);          % # unique groups
else
   cda = 0;  % CCorA:
   grp = []; % grouping vector not provided
end

if ~isempty(W)
   if size(W,2) ~= size(Y,2)
      error('W & Y must have the same # columns!');
   end
end
% ----------------------

X = f_center(X); % center on column mean
r = size(X,2);   % # of groups-1 (CDA) or # explanatory variables (CCorA)

% Hat (Projection) matrix:
% H = X*(inv(X'*X))*X';
%
% Find H without taking inverse using QR decomposition & some fancy
% algebraic substitution (courtesy news://comp.soft-sys.matlab):
[B,R1] = qr(X,0);
H      = B*B';

% Symmetric dissimilarity matrix of response variable:
yDis = f_dis(Y,dis);

% Gower's centered matrix:
A   = -0.5*(yDis.^2);
I   = eye(n,n);
uno = ones(n,1);
G   = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno'));

% Principal Coordinates:
[Q,Qevals] = f_eig(G);

% Only need n-1 axes for n objects:
Q      = Q(:,1:(end-1));
Qevals = Qevals(1:(end-1));

% Discard eigenvalues = 0:
tol    = sqrt(eps);
idx    = find (abs(Qevals)>tol);
Q      = Q(:,idx);
Qevals = Qevals(idx);

% Percent Variation Explained:
Qexpl = (Qevals / sum(Qevals)) * 100;
Qexpl = [Qexpl cumsum(Qexpl)];

% Remove PCoA's associated with neg eigenvalues and/or explains > 100%:
idx    = (Qevals>0 & Qexpl(:,2)<=100);
Q      = Q(:,idx);
Qevals = Qevals(idx);
Qexpl  = Qexpl(idx,:);
noAxes = numel(Qevals); % # retained axes

% -----Get range of m to try or ask user input:-----
if (m<0) % Get range of m for f_capOptimal:
   m        = find(Qexpl(:,2)>=55 & Qexpl(:,2)<=100);
   result.m = m;
   return;
   
elseif (m==0) % Get m interactively if not specified:
   % Show eigenvalues and '% variance explained':
   fprintf('\n-------------------------------\n');
   fprintf('        Matrix Q: \n')
   fprintf('-------------------------------\n');
   fprintf('Axis: Eigenvalues: Explained:\n')
   for i = 1:noAxes
      fprintf('%2d: %10.4f %10.2f%% %10.2f%% \n',i,Qevals(i),Qexpl(i,1),Qexpl(i,2));
   end
   fprintf('\n-------------------------------\n');
   
   % Get user to specify # of axes to retain:
   inputText = ['Specify # of Axes in Q to retain (1-' num2str(noAxes) ') ? '];
   m         = input(inputText);
end

% Avoid potential errors in LOO procedure when using Euclidean distance:
if (m>noAxes)
   m = noAxes;
end
% --------------------------------------------------

% Subset of Q:
Qm      = Q(:,1:m);
QmEvals = Qevals(1:m);

% Eigenvalue decomposition:
[U,evals]       = f_eig(Qm'*H*Qm);     % fitted Y
[Ures,evalsRes] = f_eig(Qm'*(I-H)*Qm); % residuals

% Project PCoA's in canonical space (L&L,1998):
siteScores = Qm*U;                % site scores in space of Y (= Qstar) (eq.11.12 in L&L:98)
fitScores  = (Qm'*H*Qm)*U;        % fitted site scores in space of X    (eq.11.13 in L&L:98)
resScores  = (Qm'*(I-H)*Qm)*Ures; % residual fit scores in space of residuals (fig.11.2)

% See L&L:98 page 591 for siteScores of non-canonical axes

% Scale canonical axes to sqrt of their eigenvalue:
siteScores = siteScores .* repmat(abs((evals.^0.5)'),size(siteScores,1),1);
fitScores  = fitScores  .* repmat(abs((evals.^0.5)'),size(fitScores,1),1);
resScores  = resScores  .* repmat(abs((evalsRes.^0.5)'),size(resScores,1),1);

% Keep copy of untrimmed axes for W:
if ~isempty(W)
   wU     = U;
   wEvals = evals;
end

% Trim unnecessary axes:
s          = min([r m (n-1)]); % # non-zero canonical eigenvalues
sRes       = min([m (n-1)]);   % # non-zero residual eigenvalues
U          = U(1:s,1:s);
evals      = evals(1:s);
evalsRes   = evalsRes(1:sRes);
siteScores = siteScores(:,1:s); % canonical axes
fitScores  = fitScores(:,1:s);
resScores  = resScores(:,1:sRes);

% Only keep residual axes with variances larger than 0 (L&L,1998:p.590)
idx           = find(var(resScores) < 0.000001 == 1);
evalsRes(idx) = [];
Ures(:,idx)   = [];

% Stats:
trc = sum(evals); % canonical trace statistic (= t2)
grs = evals(1);   % greatest root statistic   (= t3)

% Correlation of original Y variables with each canonical axis:
if (opt<1)
   yCor = zeros(ncY,s); % preallocate
   for i = 1:s
      for j = 1:ncY
         yCor(j,i) = f_corr(Y(:,j),siteScores(:,i));
      end
   end
else
   yCor = NaN; % skip these calculations when f_cap is called from f_capOptimal
end

% Correlation of X variables with each canonical axis:
if (cda==0) || (opt<0) % CCorA
   xCor = zeros(1,s); % preallocate
   for i = 1:s
      xCor(i) = f_corr(X,fitScores(:,i)); % SEE PAGE 591 in L&L:98
   end
else % CDA:
   xCor = NaN;
end

% Proportion of variance:
totVar = sum([evals;evalsRes]); % total variance in Y (= total inertia)
canVar = evals/totVar;          % fraction of variance explained by canonical axes
resVar = evalsRes/totVar;       % fraction of variance explained by residual axes

% Find group centroids in canonical space if CDA:
if (cda>0)
   if (sm>0)
      [centroids,nul,mDis] = f_centroid(siteScores,grp,2); % spatial median
   else
      [centroids,nul,mDis] = f_centroid(siteScores,grp,1); % centroid
   end
else
   centroids = NaN;
end
clear nul;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Classify unknowns:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(W)
   nW           = size(W,1);          % # obs in W
   wQm          = nan(nW,size(Qm,2)); % preallocate
   RSS          = nan(nW);
   siteScores_W = nan(nW,s);
   grpW         = nan(nW,1);
   PP           = nan(nW,nGrp);
   mDisW        = nan(nW,1);
   
   % Distance between Y & W:
   wDis = f_dis([Y;W],dis);
   
   % Remove inter-Y & inter-W distances:
   wA = -0.5*(wDis(n+1:end,1:n).^2); clear wDis;
   
   % Center wA according to G (eq. 7 in Anderson & Robinson, 2003
   % & eq. A.1/A.7 in Anderson & Willis, 2003):
   wG = wA - repmat(mean(A),nW,1) - repmat(mean(wA,2),1,size(wA,2))...
      + repmat(mean(A(:)),size(wA));
   
   % Classify unknown observations:
   for i = 1:nW
      % Get PCoA's of W (= wQm) (eq. A.9 in Anderson & Willis, 2003):
      % Note: wQm = (wG'*Qm)./repmat(evals',n,1) works only if wG & Qm are both
      % square matices of the same size, so wG is replicated to make it square
      wQm_temp = (repmat(wG(i,:),n,1)*Qm) ./ repmat(QmEvals',n,1);
      
      % Trim replicated duplicates, collect results:
      wQm(i,:) = wQm_temp(1,:); clear wQm_temp;
      
      % LOO Residual Sum-of-Squares (eq. A.11 in Anderson & Willis 2003):
      if ~isempty(b)
         RSS(i) = trace((b(i) - wQm(i,:)*Qm'*B)'*(b(i) - wQm(i,:)*Qm'*B));
      else
         RSS = NaN;
      end
      
      % Project wQm into canonical space (= wQstar), scale by eigenvalues:
      wQstar = wQm(i,:)*wU .* abs(wEvals.^0.5)';
      
      % Trim unnecessary axes:
      siteScores_W(i,:) = wQstar(:,1:s);
      
      % Classify each of W's siteScores to group with closest centroid:
      D = sum((repmat(siteScores_W(i,:),nGrp,1) - centroids).^2,2)'; % squared distances
      L              = 1 - (D / max(D(:))); % convert distance to likelihood
      PP(i,:)        = L / sum(L);          % posterior probabilities
      % [null,grpW(i)] = max(PP(i,:),[],2); % classify W obs
      [null,idxC]     = max(PP(i,:),[],2);  % get index to maximum PP
      grpW(i)         = idxC;               % classify W obs based on max PP
      
      % Is distance-to-nearest-centroid less than/or equal to the maximum
      % distance-to-centroid making up that group:
      mDisW(i) = (min(D) <= mDis(idxC));  
   end
   
   % Calculate margins (= confidence in each classification; after Beiman &
   % Cutler, 2003: p. 14):
   PP_S = sort(PP,2,'descend');  % sort PP's descending, by column
   marg = PP_S(:,1) - PP_S(:,2); % difference between 2 largest PP's
else
   wQm          = NaN;
   RSS          = NaN;
   siteScores_W = NaN;
   grpW         = NaN;
   PP           = NaN;
   marg         = NaN;
   mDisW        = NaN;
end
% ----------------------------


%-----Permutation tests:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   randTrc = zeros(iter-1,1); % preallocate
   randGrs = zeros(iter-1,1); % preallocate results array
   
   for i = 1:(iter-1) % observed value is considered a permutation
      
      QmPermed = f_shuffle(Qm,4); % permute obs (rows)
      
      [null,evalsPermed] = f_eig(QmPermed'*H*QmPermed); % Eigenanalysis
      randTrc(i) = sum(evalsPermed(1:s)); % permuted trace stat
      randGrs(i) = evalsPermed(1);        % permuted greatest root stat
      
   end
   j1 = find(randTrc >= trc); % get randomized stats >= to observed statistic
   j2 = find(randGrs >= grs);
   
   p.trace = (length(j1)+1)./(iter); % count values & convert to probability
   p.grs   = (length(j2)+1)./(iter);
else
   p.trace = NaN;
   p.grs   = NaN;
end
%-----------------------------


% -----Leave-One-Out Cross-validation:-----
if (loo>0)
   fprintf('\nPerforming LOO CV with m = %d...\n',m);
   
   % Predicted site scores & group membership for OMITTED obs:
   cls_o = nan(size(grp)); % preallocate
   RSS_o = nan(n,1);
   for omit = 1:n
      idx         = 1:n;                   % index of all Y's observations
      idx(omit)   = [];                    % leave one observation out
      result_o    = f_cap(Y(idx,:),dis,grp(idx),Y(omit,:),sm,0,0,m,0,B(omit,:),1); % loo=0
      cls_o(omit) = result_o.grpW;         % predicted group of omitted obs
      RSS_o(omit) = result_o.RSS;          % Residual Sum-of-Squares for omitted obs
   end
   
   % Leave-One-Out classification error rates:
   loo_err = f_errRate(grp,cls_o);
   
   % Leave-One-Out Residual Sums-of-Squares:
   loo_RSS = sum(RSS_o);
   
   % Akaike's Information Criterion for classification: <<<<< USE RSS if this is CCorA
   AIC = f_AIC(loo_err.tot,n,m,0,1);
   
else
   loo_err = NaN;
   loo_RSS = NaN;
   AIC     = NaN;
end
% -----------------------------------------



% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   if (cda>0)
      fprintf(' CAP - Canonical Discriminant Analysis:\n');
   else
      fprintf(' CAP - Canonical Correlation Analysis: \n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf('Trace Stat    = %-3.4f  p =  %3.5f \n',trc,p.trace);
   fprintf('Greatest Root = %-3.4f  p =  %3.5f \n',grs,p.grs);
   fprintf('No. of permutations = %d \n',iter);
   fprintf('--------------------------------------------------\n');
   fprintf('No. of axes of Q used (m)     = %d \n',m);
   fprintf('Variability of yDis explained = %-3.2f %s \n',Qexpl(m,2),'%');
   fprintf('Canonical Correlations:\n');
   fprintf('  %-3.4f',evals.^0.5);
   fprintf('\n')
   fprintf('Squared Canonical Correlations (= delta^2):\n');
   fprintf('  %-3.4f',evals);
   fprintf('\n==================================================\n');
end


if ((verb>0) && (loo>0) && (cda>0))
   fprintf('\n==================================================\n');
   fprintf('            LOO CROSS-VALIDATION\n'                 );
   fprintf('            Classification Success: \n'               );
   fprintf('--------------------------------------------------\n' );
   
   fprintf('Group        Correct  \n');
   for j=1:nGrp
      fprintf('%s %d %s %10.1f %s \n','  ',uGrp(j),'     ',(1-loo_err.grp(j))*100,'%');
   end
   fprintf('\n\n');
   fprintf('Total Correct  = %4.2f %s \n',(1-loo_err.tot)*100,'%');
   fprintf('Total Error    = %4.2f %s \n',loo_err.tot*100,'%');
   fprintf('\n--------------------------------------------------\n' );
   fprintf('     Confusion Matrix (%s): \n','%')
   hdr = sprintf('%-6.0f ',uGrp(:)');
   fprintf(['group: ' hdr]);
   fprintf('\n')
   for j=1:nGrp
      txt = [sprintf('%6.0f ',uGrp(j)) sprintf('%6.1f ',loo_err.conf(j,:)*100)];
      fprintf(txt);
      fprintf('\n')
   end
   fprintf('\n==================================================\n\n');
end
% ---------------------------------




%-----Wrap results up into a structure:-----
result.grp          = grp_old;
result.Qm           = Qm;    % retained PCoA's
result.Qexpl        = Qexpl; % variation of yDis explained by each PCoA

result.siteScores   = siteScores; % coordinates in canonical space (= Qstar)
result.fitScores    = fitScores;  % fitted site scores
result.resScores    = resScores;  % residal site scores
result.centroids    = centroids;  % group centroids (or spatial medians) defined in X

result.yCor         = yCor; % correlation of canonical axes with original Y variables
result.xCor         = xCor; % correlation of canonical axes with X variables

result.U            = U;        % canonical eigenvectors
result.evals        = evals;    % canonical eigenvalues (= squared canonical correlations;
%                                 1st value is greatest root statistic)
result.Ures         = Ures;     % residual eigenvectos
result.evalsRes     = evalsRes; % residual eigenvalues

result.trc          = trc; % trace statistic
result.grs          = grs; % greatest root statistic
result.p            = p;   % structure of randomized probabilities: p.trace (trace)
%                            and p.grs (greatest root)

result.canVar       = canVar; % fraction of variation explained by canonical axes
result.resVar       = resVar; % residual variance

result.wQm          = wQm;          % predicted PCoA's of W
result.RSS          = RSS;          % Residual Sum-of-Squares
result.siteScores_W = siteScores_W; % predicted site scores of W
result.grpW         = grpW;         % predicted group membership of W
result.PP           = PP;           % posterior probabilities of group prediction
result.marg         = marg;         % relative measure of certainty of predictions
result.mDisW        = mDisW;        % min distance-to-centroid doesn't exceed max distance within group (= 1)

result.loo_err      = loo_err;      % Leave-One-Out classification error rates
result.loo_RSS      = loo_RSS;      % Leave-One-Out Residual Sum-of-Squares
result.AIC          = AIC;          % Akaike's Information Criterion
if (cda>0) % central tendency
   if (sm>0)
      result.ctr = 'spatial median';
   else
      result.ctr = 'centroid';
   end
else
   result.ctr    = 'none';
end
if (cda>0)
   result.method = 'CDA';
else
   result.method = 'CCorA';
end

