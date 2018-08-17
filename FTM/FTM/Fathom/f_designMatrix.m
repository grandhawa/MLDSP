function [design,interact] = f_designMatrix(grps,trim)
% - create ANOVA design matrix using dummy variables
%
% USAGE: [design,interact] = f_designMatrix(grps,trim);
%
% grps     = matrix of integers specifying factor levels or group membership;
%            (rows = observations, cols = factors)
% trim     = omitt last column from each term to avoid a singular matrix
%            (default = 1)
%
% design   = design matrix
% interact = interaction term(s); formatted as a cell array specifying 1st
%            & 2nd order interaction terms when 'GRPS' specifies 3 factors.
% 
% SEE ALSO:  f_xMatrix, f_helmert, f_modelMatrix


% -----Notes:-----
% ANOVA design matrices are used in Multiple Linear Regression, MANOVA,
% and Canonical Discriminant Analysis to specify group membership or factor
% levels. The number of columns in a design matrix is usually 1 less than
% the number of factor levels of each term. This is all you need to adequately
% represent those levels and, more importantly, allows you to take the inverse
% of (design'*design) which would otherwise be singular. TRIM is set as the
% default for this reason.
%
% For many applications the first column of the design matrix should be a row
% of 1's specifying the intercept term. This program does NOT add an intercept.
% Note that some (most?) programs automatically insert an intercept term in
% a design matrix, so check before using the output of this program.

% -----References:-----
% Neter, J., W. Wasserman, and M. H. Kutner. 1989. Applied Linear Regression
%  Models, 2nd Edition. Irwin, Boston, MA. 667 pages.

% -----Testing:-----
% This program has been tested against SAS's PROC GLMMOD.

% -----Dependencies:-----
% requires dummyvar from the Matlab Statistics Toolbox

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 8-May-02: overhauled to support 3 factor designs
% Nov-2002: added f_recode to prevent mistakes by dummyvar
% Mar-2008: changed & to &&

if (nargin  < 2), trim = 1; end;                % trim by default
if (nargout > 1), intx = 1; else intx = 0; end; % compute interaction terms

[nr,nc] = size(grps); % nc = # of factors

if (nc<2) && (intx>0)
   error('Can''t compute interaction of only ONE term!');
end

if (nc>3)
   error('Only 3 factors are currently supported');
end

for i = 1:nc
  	dum{i} = dummyvar(f_recode(grps(:,i))); % recode as dummy variables
end;

% -----Calculate interaction terms:-----
if (nc==2) && (intx>0) % 2 Factors:
   interact = sub_interaction(dum{1},dum{2});  
elseif (nc==3) && (intx>0) % 3 Factors:
   % 1st order interactions:
   interact{1} = sub_interaction(dum{1},dum{2}); % 1x2 
   interact{2} = sub_interaction(dum{1},dum{3}); % 1x3 
   interact{3} = sub_interaction(dum{2},dum{3}); % 2x3 
   % 2nd order interaction:
   interact{4} = sub_interaction(interact{1},dum{3}); % 1x2x3
end

% add dummy variables for each factor:
design = [];
for q = 1:nc
   design = [design dum{q}];
end

% omitt last column, so not a singular matrix:
if (trim>0)
   design = design(:,1:(size(design,2)-1));
   if (intx>0) && (nc==2)
      interact = interact(:,1:(size(interact,2)-1));
   elseif (intx>0) && (nc==3)
      for i=1:4
         interact{i} = interact{i}(:,1:(size(interact{i},2)-1));
      end
   end
end

% Shouldn't the interaction terms be calculated on the 'trimmed' matrices?




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iTerm = sub_interaction(A,B)
% - calculate interaction term b/n 2 factors

aCols = size(A,2); % # columns in 1st factor
bCols = size(B,2); % # columns in 2nd factor

aCounter = 0; % initialize column counters:
bCounter = 0;
cCounter = 0;

for m = 1:aCols
   aCounter = aCounter + 1; % increment counter
   for n = 1:bCols
      bCounter = bCounter + 1; % increment column counter
      cCounter = cCounter + 1; 
      iTerm(:,cCounter) = A(:,m) .* B(:,n); % interaction term
   end
   bCounter = 0; % reset column counter
end
