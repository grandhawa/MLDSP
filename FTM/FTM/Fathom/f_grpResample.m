function [R rGrp] = f_grpResample(X,grp,n,b,c)
% - within-group resampling or bootstrapping (variable size)
%
% USAGE: [R rGrp] = f_grpResample(X,grp,n,b,c)
%
%  X   = input data (rows = observations, cols = variables)
% grp  = column vector of integers specifying group membership
% n    = number of new observations per groups
% b    = resample with replacement (= bootstrapping)               (default = 0)
% c    = when b=1, perform separately across columns               (default = 0)
%
% R    = resampled version of X
% rGrp = resampled version of grp
%
% SEE ALSO: f_grpBoot, f_boot, f_shuffle

% -----Author:-----
% by David L. Jones, Oct-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 4), b =  0; end % default no bootstraping
if (nargin < 5), c =  0; end % default don't bootstrap separately across cols

grp     = grp(:); % force col vector
[nr nc] = size(X);

if (nr ~= numel(grp))
   error('# of rows in X and GRP must be equal !');
end
% -------------------------------------

uGrp  = f_unique(grp);          % unique groups, unsorted
noGrp = size(uGrp,1);           % # of groups
R     = repmat(NaN,noGrp*n,nc); % preallocate

% Build resampled vector specifying group membership:
rGrp  = [];
for i = 1:noGrp
   rGrp = [rGrp; repmat(uGrp(i),n,1)];
end

% Resample each group separately:
for i=1:noGrp
   idxX = find(grp  == uGrp(i));
   idxR = find(rGrp == uGrp(i));
   
   if ((numel(idxX) < n) && (b==0))
      error(['Group ' num2str(i) ' is too small to resample N rows WITHOUT replacement!']);
   end
   
   if (b>0)
      idxY      = f_boot(idxX,c,n); % bootstrap indices to rows
      R(idxR,:) = X(idxY,:);        % resample rows of this groups
   else
      idxY      = f_shuffle(idxX);  % shuffle indices to rows
      idxY      = idxY(1:n);        % subsample indices
      R(idxR,:) = X(idxY,:);        % resample rows of this groups
   end
end


