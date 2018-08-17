function result = f_nnMLPcv(X,grp,hid,verb)
% - leave-one-out cross validation for f_nnMLP
%
% USAGE: result = f_nnMLPcv(X,grp,hid,verb)
%
% X     = matrix of training data (rows = obs, cols = variables) 
% grp   = vector of integers specifying group membership of X
% hid   = # hidden neurons                        (default = 5)
% verb  = verbose display of results              (default = 1);
% 
% SEE ALSO: f_nnMLP, f_nnMLP632, f_pnnCV, f_cdaCV

% -----Notes:-----
% This function is used to diagnose an MLP Feedforward Neural Network using
% the Leave-One-Out method of Cross-Validation.

% -----Author:-----
% by David L. Jones, Mar-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 3), hid   = 5;  end; % default 5 hidden neurons
if (nargin < 4), verb  = 1;  end; % default verbose output

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% ---------------------------------------

n     = size(X,1);    % # of rows (observations)
X     = f_stnd(X);    % stardardize columnwise
uGrp  = unique(grp);  % unique groups
noGrp = size(uGrp,1); % # of groups

% preallocate:
cls = repmat(NaN,n,1);

% Notify user:
if verb>0
   fprintf('Computing Leave-One-Out statistics for %d observations... \n',n);
end

for omit = 1:n
   if verb>0
      fprintf(' -> classifying observation %d \n',omit);
   end
   
   idx       = [1:n]; % index of all observations
   idx(omit) = [];    % leave one observation out
   
   [p,t,net]  = f_nnMLP(X(idx,:),grp(idx),hid);     % create net   
   [net_2,tr] = train(net,p,t);                     % train  net
   cls(omit)  = f_compet([sim(net_2,X(omit,:)')]'); % classify omitted observation
   
   % clean up:
   clear p t net net_2 tr;
      
end

% Classification error rates:
err = f_errRate(grp,cls);


% -----Send output to display:-----
if (verb>0)
   
   fprintf('\n==================================================\n');
   fprintf('            f_nnMLP CROSS-VALIDATION\n'                 );
   fprintf('            Classification Success: \n'               );
   fprintf('--------------------------------------------------\n' );
   
   fprintf('Group        Corrrect  \n');
   for j=1:noGrp
      fprintf('%s %d %s %-10.1f %s \n',['  '],uGrp(j),['     '],[1-err.grp(j)]*100,['%']);
   end
      
   fprintf('\n\n');
   fprintf('Total Correct  = %-3.1f %s \n',(1-err.tot)*100,['%']);
   fprintf('Total Error    = %-3.1f %s \n',err.tot*100,['%']);
   
   fprintf('\n--------------------------------------------------\n' );
   fprintf('     Confusion Matrix (%s): \n',['%'])
      hdr = sprintf('%-6.0f ',uGrp(:)');
   fprintf(['group: ' hdr]);
   fprintf('\n')
   for j=1:noGrp
      txt = [sprintf('%-6.0f ',uGrp(j)) sprintf('%-6.1f ',err.conf(j,:)*100)];
      fprintf(txt);
      fprintf('\n')
   end
   fprintf('\n==================================================\n\n');
end


