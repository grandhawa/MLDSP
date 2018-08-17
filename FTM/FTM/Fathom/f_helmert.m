function D = f_helmert(x)
% - Helmert orthogonal contrast codes for an ANOVA design matrix
% 
% USAGE: D = f_helmert(x)
% 
% x = column vector of integers specifying factor levels or group membership
% D = design matrix
% 
% SEE ALSO: f_designMatix, f_modelMatrix

% -----Notes:-----
% This method of coding requires a BALANCED design. Orthogonal coding for the
% interaction between two factors can be obtained by multiplying the Helmerth
% codes for those factors together (see exampleDesign.m).

% -----References:-----
% Legendre, P. and M. J. Anderson. 1999. Distance-based redundancy analysis:
%   testing multispecies responses in multifactorial ecological experiments.
%   Ecological Monographs 69(1): 1-24. (See Appendix C)

% -----Author:-----
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (size(x,2) > 1)
   error('X must specify a single variable & have only 1 column!')
end

x  = f_recode(x); % recode as consecutive integers
r  = size(x,1);   % # rows 
Rx = unique(x);   % treatment levels
n  = size(Rx,1);  % number of treament levels 

% Check for a balanced design:
nReps = zeros(n,1);
if max(nReps) > min(nReps)
   error('This method requires a BALANCED design!')
end

% Build a generalized Helmert contrast for n treatment levels:
H       = eye(n,n-1);         % initialize
H(H==1) = n-1:-1:1;           % populate the diagonal
H(~triu(ones(size(H)))) = -1; % replace lower tridiagonal with -1's

% Construct matrix of orthogonal contrast codes:
D     = zeros(r,n-1); % preallocate
for i=1:n
   idx      = find(x==Rx(i));
   D(idx,:) = repmat(H(i,:),size(idx,1),1);
   nReps(i) = size(idx,1);
end

