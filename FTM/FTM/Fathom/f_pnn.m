function result = f_pnn(X,grp,Y,prior,sm)
% - Probabilistic Neural Network
%
% USAGE: result = f_pnn(X,grp,Y,prior,sm)
%
% X     = matrix of training data (rows = obs, cols = variables) 
% grp   = vector of integers specifying group membership of X
% Y     = matrix of test data to predict group membership
% prior = prior probabilities; 'empirical' (default) or 'equal'
% sm    = smoothing factor of radial basis functions   (default = 0.1)
% verb  = optionally display results                   (default = 1)
%
% result        = structure having the following fields:
% result.pred   = predicted group membership
% result.Bayes  = Bayes posterior probabilities
% result.target = target probabilities (adjusted by prior)
% result.prior  = prior probabilities used to adjust Bayes

% -----Notes:-----
% This function is used to create a Probabilistic Neural Network in order
% to perform Kernel Discriminant Analysis. This is a neural net
% implementatin of nonparametric, nonlinear discriminant analysis which
% does not require the data to be multivariate normal, have homogeneous
% variances among groups, or equal group size.
%
% A SOFTMAX transfer function is used to output Bayes posterior
% probabilities. Prior probabilities (PRIOR) are used to adjust the Bayes
% probabilities according to group size when prior = 'empirical'; this
% improves the error rates for the rarer groups. No adjustments are made
% when prior = 'equal'. 

% -----References:-----
% Matlab NNET Toolbox reference &
% news://comp.ai.neural-nets

% -----Dependencies:-----
% This program requires the Matlab NNET Toolbox.

% -----Author:-----
% by Dave Jones, Jan-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input :-----
if (nargin<4), prior = 'empirical'; end; % default priors
if (nargin<5), sm    = 0.1;         end; % default spread

if (size(X,2) ~= size(Y,2))
   error('X & Y must have same # of columns!');
end

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% -------------------------------------

[rX,cX] = size(X);      % # obs & var in X
[rY,cY] = size(Y);      % # obs & var in Y
uGrp    = unique(grp);  % unique groups
noGrp   = size(uGrp,1); % # groups

% Stardardize to equally weight variables:
Y  = f_stnd(Y,X); % stardardize Y by X
X  = f_stnd(X);

% Create a new PNN:
net = newpnn(X',ind2vec(grp'),sm);

% Change transfer function:
net.layers{2}.transferFcn = 'softmax'; % Bayes posterior probabilities

% Run simulation:
Bayes = sim(net,Y');

% -----Adjust Bayes by prior probabilities:-----
switch prior     
   case 'equal'
      % no adjustment is made
   
   case 'empirical'    
      pprob = zeros(1,noGrp);
      for i=1:noGrp
         idx      = find(grp==uGrp(i));
         grpSize  = size(idx,1);
         pprob(i) = grpSize/rX;
      end
      pprob = [repmat(pprob,rY,1)]';
      Bayes = [Bayes ./ pprob]';              % divide by ratio of grp size
      Bayes = [Bayes./repmat([sum(Bayes')'],1,size(Bayes,2))]'; % renormalize
            
   otherwise
      error('Prior probabilities must be ''empirical'' or ''equal''!');
end
% ---------------------------------------------

% Predict group membership:
pred = vec2ind(compet(Bayes));

% Target (= expected probabilities of training data):
target = zeros(rY,noGrp);
for i=1:noGrp
   idx = find(grp==uGrp(i));
   target(idx,i) = 1;
end

% -----Wrap results up into a structure:-----
result.pred   = pred';
result.Bayes  = Bayes';
result.target = target;

switch prior
   case 'equal'
   result.prior = prior;
otherwise
   result.prior = [pprob(:,1)]';
end
% -------------------------------------------

