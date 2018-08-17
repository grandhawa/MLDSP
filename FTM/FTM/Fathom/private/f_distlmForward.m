function [fx,f_labels] = f_distlmForward(x,n,I,G,labels,iter)
% - utility function called by f_distlm

% -----Input/Output:-----
% x = matrix of explanatory variables
% n = # rows/colums in distance matrix
% I = I matrix
% G = Gower's centered matrix
% labels = variable labels
%
% fx       = matrix of predictor variables, in decreasing order of
%            contribution to Total Sum-of-Squares 
% f_labels = corresponding variable labels (re-ordered)

% -----Notes:-----
% This function performs a forward selection of explanatory variables using
% the criterion of maximum proportion of Sum-of-Squares explained.

% -----Author:-----
% by David L. Jones, Aug-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


% ========================================
%       Forward Variable Selection:
% ========================================

noVars = size(x,2); % # of explanatory variables
SST    = trace(G);  % Sum-of-Squares Total (McArdle & Anderson, 2001)

fx       = [ones(n,1)];           % initialize with intercept term
f_labels = {'null',size(labels)}; % preallocate

for j = 1:noVars
   
   ncX = size(x,2); % # of (remaining) explanatory variables
   for i = 1:ncX    % select a variable that maximizes R2
      zz      = [fx x(:,i)];         % add variables sequentially
      [Q1,R1] = qr(zz,0); H= Q1*Q1'; % Hat-matrix
      SSR     = trace(H*G*H);        % SS Regression (McArdle & Anderson, 2001)
      R2(i)   = SSR/SST;             % proportion of SST explained by SSR (Neter et al.,1996 [eq. 6.40])
      
   end
   
   [null,idx]  = max(R2');      % get index of variable that maximizes R2 
   fx          = [fx x(:,idx)]; % permanently add that variable to regression
   f_labels(j) = labels(idx);   % maintain corresponding variable labels
   
   % clear for next iteration:
   x(:,idx)    = []; % remove variable from further selection
   labels(idx) = []; % ditto
   clear prop;
   
end

fx(:,1)  = [];          % remove intercept term
f_labels = f_labels(:); % force column vector

