function res = f_pnnCV(X,grp,sm,verb)
% - leave-one-out cross validation for PNN
%
% USAGE: f_pnnCV(X,grp,sm,verb);
%
% X      = matrix of training data (rows = obs, cols = variables) 
% grp    = vector of integers specifying group membership of X
% sm     = smoothing factor of radial basis functions   (default = 0.1)
% verb   = optionally display results                   (default = 1)
%
% res    = total classification success rate

% -----Notes:-----
% This program is used to perform leave-one-out cross validation for a
% Probabilistic Neural Network (PNN).

% -----Author:-----
% by David L. Jones, Jan-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin<3), sm   = 0.1; end; % default smoothing factor
if (nargin<4), verb = 1;   end; % default show output

grp = grp(:); % force column vector
if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% -------------------------------------

% Get size of input:
n     = size(X,1);    % # of rows (observations)
uGrp  = unique(grp);  % unique groups
noGrp = size(uGrp,1); % # of groups

% Preallocate:
classifiedGrp = repmat(NaN,n,1);
classRateGrp  = repmat(NaN,noGrp,1);

% Notify user:
if verb>0
   fprintf(' -> Computing Leave-One-Out statistics for %d observations... \n',n);
end

% Perform Leave-One-Out:
for omit = 1:n
   idx       = [1:n]; % index of all observations
   idx(omit) = [];    % leave one observation out
   temp      = f_pnn(X(idx,:),grp(idx),X(omit,:),'empirical',sm); % PNN
   classifiedGrp(omit) = temp.pred;                               % collect results
end

% -----Classification success:-----
correct = logical([grp - classifiedGrp] ==0); % 1=pass, 0=fail

% Total classification success:
classRate = [sum(correct)/size(correct,1)]*100;

% Classification success by group:
for i=1:noGrp
   idx = find(grp==uGrp(i));
   classRateGrp(i) = [sum(correct(idx)/size(correct(idx),1))]*100;
end
% ------------------------------

res = classRate/100;


% -----Send output to display:-----
if (verb>0)
   
   fprintf('\n==================================================\n');
   fprintf('            PNN CROSS-VALIDATION\n'                 );
   fprintf('            Classification Success: \n'               );
   fprintf('--------------------------------------------------\n' );
   
   fprintf('Group        Corrrect  \n');
   for j=1:noGrp
      fprintf('%s %d %s %-10.4f %s \n',['  '],uGrp(j),['     '],classRateGrp(j),['%']);
   end
      
   fprintf('\n\n');
   fprintf('Total Correct  = %-3.4f %s \n',classRate,['%']);
   fprintf('Total Error    = %-3.4f %s \n',[100 - classRate],['%']);
   fprintf('Prior prob     = group size \n');
   fprintf('\n==================================================\n\n');
   
end
