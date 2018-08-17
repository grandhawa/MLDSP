function sol = f_pnn632(X,grp,sm,iter,verb)
% - estimate PNN error using bootstrap .632+ rule
%
% USAGE: sol = f_pnn632(X,grp,sm,iter,verb)
%
% X     = matrix of training data (rows = obs, cols = variables) 
% grp   = vector of integers specifying group membership of X
% sm    = smoothing factor of radial basis functions   (default = 0.1)
% iter  = # iterations for bootstrap resampling        (default = 50)
% verb  = verbose display of results                   (default = 1);
% 
% sol   = structure of results with the following fields:
% sol.err         = apparent error
% sol.err_1       = bootstrap error
% sol.err_632     = .632 bootstrap error
% sol.err_632plus = .632+ bootstrap error
% sol.R           = relative overfitting rate;
% sol.gamma       = no-information rate;
% sol.p           = prior probabilities;
% sol.q           = posterior probabilities;
%
% SEE ALSO: f_pnn, f_pnnSm, f_pnnCV, f_nnMLP632, f_cda632

% -----References:-----
% Efron, B. and R. Tibshirani. 1997. Improvements on cross-validation: the
%   .632+ bootstrap method.
% Furlanello, C., S. Merler, C. Chemini, and A. Rizzoli. An application of
%   the bootstrap 632+ rule to ecological data.
% Ambroise, C. and G. J. McLachlan. 2002. Selection bias in gene extraction 
%   on the basis of microarray gene-expression data. Proc. Natl. Acad. Sci.
%   USA 99(10): 6562-6566.
%
% Equations refer to Efron & Tibshirani (1997) unless indicated otherwise.

% -----Author:-----
% by David L. Jones, Feb-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 3), sm    = 0.1; end; % default 50 bootstraps
if (nargin < 3), iter  = 50;  end; % default 50 bootstraps
if (nargin < 4), verb  = 1;   end; % default verbose output

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% ---------------------------------------

n = size(X,1); % # of observations
X = f_stnd(X); % stardardize columnwise

% Create, train, & classify PNN:
temp = f_pnn(X,grp,X,'empirical',sm);
cls  = temp.pred;

% Apparent error (training set == test set):
err = [sum(logical([cls - grp] ~= 0))]/n;

% -----No-information rate:-----
uGrp  = unique(grp);  % unique groups
noGrp = size(uGrp,1); % # groups

% preallocate:
p     = zeros(1,noGrp);
q     = p;
gamma = 0;

for i=1:noGrp % after Ambroise & McLachlan, 2002 (eq.3):
   p_idx      = find(grp==uGrp(i));    % input
   q_idx      = find(cls==uGrp(i));    % output
   p_grpSize  = size(p_idx,1);
   q_grpSize  = size(q_idx,1);
   p(i)       = p_grpSize/n;           % priors
   q(i)       = q_grpSize/n;           % posteriors
   gamma      = gamma + p(i)*[1-q(i)]; % no-information rate
end
% ------------------------------

% -----Bootstrap:-----
err_boot = zeros(iter,1); % preallocate
for b = 1:iter      
  
   % Bootstrap training set:
   [idxB,idxT] = f_boot(1:n);
      
   % Create, train, classify bootstrap PNN:
   temp = f_pnn(X(idxB,:),grp(idxB),X(idxT,:),'empirical',sm);
   cls  = temp.pred;
   
   % Mis-classification Error (training ~= test set);
   err_boot(b) = [sum(logical([cls - grp(idxT)] ~= 0))]/length(idxT); 
end   

% Bootstrap error:
err_1 = mean(err_boot(:));
% --------------------   

% Bootstrap .632 (eq. 24):
err_632 = (0.368*err) + (0.632*err_1);

% Relative overfitting rate (eq. 28):
R = (err_1 - err)/(gamma - err);

% -----Make sure R ranges from 0-1:-----
err_1 = min([err_1 gamma]);

if (err_1 > err) & (gamma > err) % (eq. 31)
   % R = R;
else
   R = 0;
end
% --------------------------------------

% Bootstrap .632+ (eq. 32):
err_632plus = err_632 + (err_1 - err)*[(0.368*0.632*R)/(1-0.368*R)];

% Wrap results up into a structure:
sol.err         = err;
sol.err_1       = err_1;
sol.err_632     = err_632;
sol.err_632plus = err_632plus;
sol.R           = R;
sol.gamma       = gamma;
sol.p           = p;
sol.q           = q;

