function result = f_ncap(Y,dis,X,grad,type,verb,m,iter,bs,opt)
% - Nonlinear Canonical Analysis of Principal Coordinates (NCAP)
%
% USAGE: result = f_ncap(Y,dis,X,'grad','type',verb,m,iter,bs);
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
% type = type of statistic to use: 'R2' (= default) or 'RDA'
% verb = optionally send results to display                        (default = 1)
% m    = # of PCoA axes to retain; m = 0 retains all               (default = 0)
% iter = # iterations for permutation test                         (default = 0)
% bs   = # iterations for bootstrapped 95% CI's                    (default = 0)
%
% result = structure of outputs with the following fields:
%  .Q      = PCoA's from dissimilarity matrix of Y
%  .Qevals = eigenvalues of Q
%  .Qexpl  = percent variation explained
%  .Qm     = subset of Q
%  .m      = # PCoA axes retained
%  .b      = optimal nonlinear canonical coefficients
%  .stat   = squared correlation coefficients (R^2) or RDA statistic
%  .type   = statistic used for the fitting criterion ('R2' or 'RDA')
%  .AIC    = Akaike's Information Criterion
%  .lin_R2 = linear CCorA squared correlation coefficients
%  .lin_b  = linear CCorA canonical coefficients
%
% SEE ALSO: f_ncapOptimal, f_cap, f_CCorA

% -----Notes:-----
% Using type = 'RDA' performs a nonlinear least squares regression fit that
% weights the PCoA's according to their relative contribution to the
% total variation in Y; using type = 'R2' uses equal weights.
%
% If the explanatory variable (= the gradient) is standardized, this must
% be taken into account in the interpretation of the results (See Lengendre
% & Legendre, 1998: p.615).

% -----References:-----
% Portions of this function are based on ideas implemented in Millar's
% 'NCAP.R' code for R.
%
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. pages 612-616
% Millar, R. B., M. J. Anderson, & G. Zunun. 2005. Fitting nonlinear
%  environmental gradients to community data: a general distance-based approach.
%  Ecology 86(8): 2245-2251.
%
% See also references in F_CAP
%
% For LSQNONLIN see: http://www.mathworks.com/support/tech-notes/1500/1504.html
% For weighted least squares fit see:
% http://www.mathworks.com/support/solutions/en/data/1-18DGY/?solution=1-18DGY

% -----Author:-----
% by David L. Jones, Mar-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% May-2014: updated documentation

% -----Set defaults & check input:-----
if (nargin < 4),  grad  = 'von'; end % default gradient if von von Bertalannfy
if (nargin < 5),  type  = 'R2';  end % default type is R^2
if (nargin < 6),  verb  = 1;     end % default send results to display
if (nargin < 7),  m     = 0;     end % default get m interactively
if (nargin < 8),  iter  = 0;     end % default 0 iterations
if (nargin < 9),  bs    = 0;     end % default no bootstrap iterations
if (nargin < 10), opt   = 0;     end % internal flag indicating call by f_ncapOptimum

% Check that gradient is univariate:
if size(X,2)>1, error('Only 1 column of X is currently supported!'); end

% Check if Y,X are compatible:
if size(Y,1) ~= size(X,1), error('Y & X need same # of rows'); end

% Check gradient:
if ~ischar(grad)
   error('GRAD should be a character array!');
end
if ~(ismember(grad,{'von' 'hyp' 'log' 'lin'}))
   error('GRAD should be ''von'', ''hyp'', ''log'', or ''lin''!');
end

% Check type:
if ~ischar(type)
   error('TYPE should be a character array!');
end
type = upper(type); % force upper case
if ~(ismember(type,{'R2' 'RDA'}))
   error('TYPE should be ''R2'' or ''RDA''!');
end

% Check permutation/bootstrap iterations
if (iter>0) && (bs>0)
   error('Set ITER=0 when generating bootstrapped 95% CI''s!');
end
% -------------------------------------

n = size(X,1); % get # of observations

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
Qexpl = (Qevals / sum(Qevals));
Qexpl = [Qexpl cumsum(Qexpl)]*100;

% Remove PCoords associated with neg eigenvalues and/or explains > 100%:
idx    = (Qevals>0 & Qexpl(:,2)<=100);
Q      = Q(:,idx);
Qevals = Qevals(idx);
Qexpl  = Qexpl(idx,:);
noAxes = numel(Qevals); % # retained axes

if (opt>0)
   % Stop here if called by f_ncapOptimum:
   if (m==0)
      result.m  = noAxes;
   else
      result.m  = m;
   end
   result.propY = Qexpl(:,2);
   return;
   % Retain all axes:
elseif (m==0)
   m = noAxes;
end

% Subset Q:
Qm = Q(:,1:m);

% Weight the PCoA's according to their relative importance:
if isequal(upper(type),'RDA')
   wts = Qexpl(1:m,1)/100;
else
   wts = ones(m,1); % equal weights
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      LINEAR CCorA:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lin    = f_CCorA(repmat(wts',n,1).*Qm,X,0,0);
lin_R2 = lin.evals; % linear squared correlation coefficients
lin_b  = lin.Cx;    % linear canonical coefficients to initialize b0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    NONLINEAR CCorA:                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use nonlinear optimization to estimate coefficients (i.e., b):
opt = optimset('MaxIter',500,'DiffMinChange',1e-6,'TolFun',1e-8,'TolX',1e-8,'Display','off');
b0  = lin_b;
LB  = zeros(size(b0)); % constrain b > 0
UB  = [];

% Optionally display the value of m used:
if (verb>0)
   fprintf('\n-----Finding nonlinear coefficients for m = %d:-----\n',m);
end

% Nonlinear optimization:
b = lsqnonlin(@obj_NLCCor,b0,LB,UB,opt); % final coefficients
   function res = obj_NLCCor(b0)
      % This is a nested objective function for LSQNONLIN; since LSQNONLIN ONLY
      % allows its objective functions to have 3 inputs, we use a 'nested'
      % version here to access the other variables we need in the parent
      % function's workspace. Note that X & Q could be passed with b to
      % LSQNONLIN (e.g., b = lsqnonlin(@obj_NLCCor,b0,LB,UB,opt,X,Q), but X
      % & Q would then be updated along with b during each iteration of the
      % fitting process; however, we instead want to optimize b given a
      % constant X & Q.
      
      % Get fitted values of Y using specified gradient:
      switch grad
         case 'von'
            Yfit = (1 - exp(-X*b0)); % (eq. 5 in Millar et al. 2005)
         case 'hyp'
            Yfit = (X*b0)/(1+X*0);
         case 'log'
            Yfit = exp(X*b0)/(1+exp(X*b0));
         case 'lin'
            Yfit = 1+(X*b0);
         otherwise
            error('Unknown GRAD!')
      end
      
      % Center data:
      Yfit = f_center(Yfit);
      
      % Sum-of-Squares of Covariance matrices (see NLCCor in Millar's 'NCAP.R'...
      % and L&L,1998:p613):
      SSyy = trace(Yfit'*Yfit);
      SSyQ = sum(wts.*diag((Yfit'*Qm)'*Yfit'*Qm));
      
      % Correlation:
      stat  = SSyQ/SSyy;
      
      % Return residuals to LSNONLIN so it can minimize their Sum-of-Squares:
      res = 1 - stat;
      
      % Monitor the fitting process:
      if (verb>0)
         fprintf('stat = %11.9f   SSres = %18.9f    b = %11.9f\n',stat,trace(res'*res),b0);
      end
   end
% --------------------------------------------------------------

% F-statistic for nonlinear vs. linear gradient:
if isequal(type,'R2')
   F = (stat - lin_R2)/(1 - stat);
else
   F = NaN;
end

% Akaike's Information Criterion:
AIC = f_AIC(1-stat,n,m,0);

% Optionally send output to display:
if (verb>0)
   fprintf('\nFirst %d axes of Q explains %3.2f%% of Y\n',noAxes,Qexpl(end,2));
end

%-----Permutation tests:-----
if (iter>0)
   if ~isequal(type,'R2')
      error('Permutation tests require TYPE = ''R2''!')
   end
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   permStat = zeros(iter-1,1); % preallocate
   permF    = zeros(iter-1,1);
   for i = 1:(iter-1) % observed value is considered a permutation
      permX       = f_shuffle(X,4); % permute obs (rows)
      permRes     = f_ncap(Y,dis,permX,grad,type,0,m,0,0,0);
      permStat(i) = permRes.stat;   % permuted R2
      permF(i)    = permRes.F;      % permuted F-statistic
   end
   % Get randomized stats >= to observed statistic
   p.stat = (sum(permStat >= stat)+1)/iter; % count values & convert to probability
   p.F    = (sum(permF    >= F)+1)/iter;
else
   p.stat = NaN;
   p.F    = NaN;
end
%-----------------------------

% -----Bootstrapped 95% CI of b:-----
if (bs>0)
   fprintf('\nBootstrapping the data %d times...\n',bs);
   conf = 0.95;
   bBoot = nan(bs,numel(b)); % preallocate
   for i=1:bs
      idxB       = randi(n,n,1); % get bootstrapped index to rows of Y,X
      boot       = f_ncap(Y(idxB,:),dis,X(idxB,:),grad,type,0,m,0,0,0);
      bBoot(i,:) = boot.b;
      clear boot;
   end
   cv = prctile(bBoot,(1-(1 - conf))*100);
   CI = [abs(b-cv) b+cv];
else
   CI = NaN;
end
% -----------------------------------

% Wrap results up into a structure:
result.Q      = Q;      % PCoA's from dissimilarity matrix of Y
result.Qevals = Qevals; % eigenvalues of Q
result.Qexpl  = Qexpl ; % percent variation explained
result.Qm     = Qm;     % subset of Q
result.m      = m;      % # PCoA axes retained
result.b      = b;      % nonlinear canonical coefficients
result.stat   = stat;   % nonlinear squared correlation coefficients (R^2) or RDA
result.type   = type;   % statistic use for the fitting criterion (R^2 or RDA)
result.AIC    = AIC;    % Akaike's Information Criterion
result.lin_R2 = lin_R2; % linear CCorA squared correlation coefficients
result.lin_b  = lin_b ; % linear CCorA canonical coefficients
result.F      = F;      % F-statistic for nonlinear vs. linear gradient
result.p      = p;      % corresponding permutation-based p-values
result.CI     = CI;     % bootstrapped 95% confidence intervals

end % final termination required for nested functions

