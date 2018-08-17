function err = f_errRate(grp,cls)
% - error rate for a classifier
%
% USAGE: err = f_errRate(grp,cls)
%
% grp = input vector of integers specifying group membership
% cls = group membership predicted by the classifier
%
% err = structure of results having the following fields:
% err.tot  = total error
% err.grp  = total error for each group
% err.uGrp = list of groups
% err.conf = confusion matrix (proportion of CLS classified as GRP)
%
% SEE ALSO: f_cdaCV, f_boot632, f_chanceClass

% -----Notes:-----
% This function is used to determine error rates for classifiers; in
% particular, it is called by f_cdaCV. More robust measures of
% generalization error may be achieved using f_boot632, etc.

% -----Author:-----
% by David L. Jones, Feb-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2011: replaced 'unique' with 'f_unique'; added uGrp to output

% -----Check input:-----
grp = grp(:);
cls = cls(:);

if size(grp,1) ~= size(cls,1)
   error('GRP and CLS must be same size!')
end
% -----------------------

n     = size(grp,1);   % # obs
uGrp  = f_unique(grp); % unique groups, unsorted
noGrp = size(uGrp,1);  % # of groups

% Preallocate:
err.tot  = NaN;
err.grp  = zeros(noGrp,1);
err.uGpr = uGrp;
err.conf = zeros(noGrp,noGrp);

for i=1:noGrp
   idx     = find(grp==uGrp(i));
   grpSize = size(idx,1);
   
   % Confusion matrix:
   for j=1:noGrp
      err.conf(i,j) = [sum(logical(cls(idx) == uGrp(j)))]/grpSize;
   end   
end

% Total error rate:
err.tot = [sum(logical([grp - cls] ~= 0))]/n; % 1=pass, 0=fail

% Error rate by group:
err.grp = diag([1 - err.conf]);
