function sol = f_cda632(X,grp,method,iter,verb)
% - estimate CDA error rate using bootstrap .632+ rule
%
% USAGE: sol = f_cda632(X,grp,method,iter,verb)
%
% X      = matrix of training data (rows = obs, cols = variables) 
% grp    = vector of integers specifying group membership of X
% method = method of classification:
%          1 : linear     (default)
%          2 : quadratic
%          3 : mahalanobis
%
% iter  = # iterations for bootstrap resampling   (default = 50)
% verb  = verbose display of results              (default = 1);
% 
% sol   = structure of results with the following fields
%  .err         = apparent error
%  .err_1       = bootstrap error
%  .err_632     = .632 bootstrap error
%  .err_632plus = .632+ bootstrap error
%  .R           = relative overfitting rate;
%  .gamma       = no-information rate;
%  .p           = prior probabilities;
%  .q           = posterior probabilities;
%
% SEE ALSO: f_cdaCV, f_cdaBCV, f_cda, classify

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
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), method = 1;  end  % linear method by default
if (nargin < 4), iter   = 50; end  % default 50 bootstraps
if (nargin < 5), verb   = 1;  end  % send output to display by default

switch method
   case 1, type = 'linear';
   case 2, type = 'quadratic';
   case 3, type = 'mahalanobis';
   otherwise, error('Classification method must be 1, 2, or 3!');
end

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% ---------------------------------------

n = size(X,1); % # of observations
X = f_stnd(X); % stardardize columnwise

% Apparent Error (training set == test set):
cls = classify(X,X,grp,type,'empirical');
err = [sum(logical([cls - grp] ~= 0))]/n;

% -----No-information rate:-----
uGrp  = unique(grp);  % unique groups
noGrp = size(uGrp,1); % # groups

% Preallocate:
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

fprintf('\nTraining %d bootstrap samples...\n',iter);
% -----Bootstrap:-----
err_boot = zeros(iter,1); % preallocate
for b = 1:iter      
   if (verb==1),fprintf(' -> %d of %d \n',b,iter); end;
   
   % Bootstrap training set:
   [idxB,idxT] = f_boot(1:n);
      
   % Classify test set (Training and Test sets don't overlap):
   cls = classify(X(idxT,:),X(idxB,:),grp(idxB),type,'empirical');
   
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

if (err_1 > err) && (gamma > err) % (eq. 31)
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

