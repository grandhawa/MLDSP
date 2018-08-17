function  [X,int] = f_xMatrix(Y,orth,nest)
% - design matrix of contrast codes for ANOVA linear models
%
% USAGE: [X,int] = f_xMatrix(Y,orth,nest)
%
% Y    = matrix of integers specifying factor levels or group membership;
%        (rows = observations, cols = factors)
% orth = create orthogonal contrast codes                       (default = 0)
% nest = factors are nested (i.e., factor 2 within factor 1)    (default = 0)
%        Eaxmple:
%        factor 1 = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]'; % main factor
%        factor 2 = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]'; % nested factor
%
% X    =  design matrix (all crossed factors combined, or 1 nested factor)
% int  = interaction term(s) for un-nested factors
% 
% SEE ALSO: f_dummy, f_helmert, f_designMatrix, f_modelMatrix

% -----References:-----
% Anderson, M. J. 2003. XMATRIX: a FORTRAN computer program for calculating
%   design matrices for terms in ANOVA designs in a linear model. Department of
%   Statistics, University of Auckland, New Zealand.

% -----Author:-----
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
[nr,nf]  = size(Y); % # obs & factors

if (nargin < 2), orth = 0; end % default to non-orthogonal codes
if (nargin < 3), nest = 0; end % default no nested factors

if (nf < 2) && (nest==1)
   error('Nested designs require at last 2 factors!');
end

if (nf > 3) && (nest==0)
   error('Only 3 crossed factors are currently supported!');
end
% -------------------------------------

% -----Main Effects:-----
if (nest>0) % Nested Factors:
   D = sub_contrast(Y(:,nf),orth); % initialize with lowest factor in hierarchy
   for i = nf:-1:2                 % move up the hierarchy
      D = sub_nested(f_dummy(Y(:,i-1),0),D);
   end
else        % Crossed Factors:
   D{nf} = []; % preallocate
   for i = 1:nf
      D{i} = sub_contrast(Y(:,i),orth);
   end
end

% -----Interaction terms:-----
if (nf==1) || (nest==1)        % 1 Factor (or Nested):
   int = NaN;

elseif (nf==2) && (nest==0)    % 2 Factors:
   int = sub_inter(D{1},D{2});

elseif (nf==3)                 % 3 Factors:
   % 1st order interactions:
   int{1} = sub_inter(D{1},D{2});   % 1x2
   int{2} = sub_inter(D{1},D{3});   % 1x3
   int{3} = sub_inter(D{2},D{3});   % 2x3

   % 2nd order interaction:
   int{4} = sub_inter(int{1},D{3}); % 1x2x3
end


%-----Combine factors:-----
if (nest>0) % Nested Factors:
   X = D;
else        % Crossed Factors:
   X = [];
   for q = 1:nf
      X = [X D{q}];
   end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = sub_contrast(Y,orth)
% - convert qualitative variables to contrast codes using [1 0 -1]
Y  = f_recode(Y); % force as consecutive integers
nr = size(Y,1);   % # rows
G  = unique(Y);   % groups
ng = size(G,1);   % # unique groups

X = zeros(nr,ng-1); % preallocate
for i = 1:ng-1
   X(Y == G(i),i)   =  1; % obs belongs to this group
   X(Y == G(i)+1,i) = -1; % obs belongs to next group
end

if (orth>0) % Make orthogonal:
   X = f_pca(X,0,1,sqrt(eps));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iTerm = sub_inter(A,B)
% - calculate interaction term b/n 2 factors
ncA    = size(A,2);
ncB    = size(B,2);
iTerm  = repmat(A,1,ncB) .* repmat(B,1,ncA);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = sub_nested(A,B)
% - B is nested within A
c = size(A,2);
N = [];
for i = 1:c
   N = [N  [repmat(A(:,i),1,size(B,2)) .* B]];
end
