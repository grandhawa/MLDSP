function sol = f_nnMLP632(X,grp,iter,hid,verb)
% - estimate MLP error rate using bootstrap .632+ rule
%
% USAGE: sol = f_nnMLP632(X,grp,iter,hid,verb)
%
% X     = matrix of training data (rows = obs, cols = variables) 
% grp   = vector of integers specifying group membership of X
% iter  = # iterations for bootstrap resampling   (default = 50)
% hid   = # hidden neurons                        (default = 5)
% verb  = verbose display of results              (default = 1);
% 
% sol   = structure of results with the following fields
% sol.err         = apparent error
% sol.err_1       = bootstrap error
% sol.err_632     = .632 bootstrap error
% sol.err_632plus = .632+ bootstrap error
% sol.R           = relative overfitting rate;
% sol.gamma       = no-information rate;
% sol.p           = prior probabilities;
% sol.q           = posterior probabilities;
%
% SEE ALSO: f_nnMLP, f_nnMLPcv, f_pnn632, f_cda632

% -----Notes:-----
% Efron & Tibshirani (1997) suggest that only 50 bootstrap samples are
% needed.
%
% APPARENT ERROR (= resubstitution error) uses the same data to train the
% net and to test it.
%
% NO-INFORMATION RATE = expected error rate of the classifier when the data
% provides no information for classification (i.e., input and output are
% independent).

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
% by Dave Jones, Feb-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 3), iter  = 50; end; % default 50 bootstraps
if (nargin < 4), hid   = 5;  end; % default 5 hidden neurons
if (nargin < 5), verb  = 1;  end; % default verbose output

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% ---------------------------------------

n = size(X,1); % # of observations
X = f_stnd(X); % stardardize columnwise

% Create & train the net:
[p,t,net]  = f_nnMLP(X,grp,hid);        % create net   
[net_2,tr] = train(net,p,t);            % train  net
cls        = f_compet([sim(net_2,p)]'); % classify training set

% Apparent Error (training set == test set):
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

fprintf('\nTraining %d bootstrap samples...\n',iter);
% -----Bootstrap:-----
err_boot = zeros(iter,1); % preallocate
for b = 1:iter      
   if (verb==1),fprintf(' -> %d of %d \n',b,iter); end;
   
   [idxB,idxT] = f_boot(1:n);                      % bootstrap training set
   [B,t,net]   = f_nnMLP(X(idxB,:),grp(idxB),hid); % create bootstrap net   
   [net,tr]    = train(net,B,t);                   % train bootstrap net

   cls         = f_compet([sim(net,X(idxT,:)')]'); % classify test set
   
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

