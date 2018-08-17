function [nPerms,pd] = f_multicomb(grps,trim,max)
% - generate permutation distribution of grouping labels
%
% USAGE: [nPerms,pd] = f_multicomb(grps,{trim},{max});
%
% grps = vector of integers specifying group membership
%        e.g., grps = [1 1 1 2 2 2 2 3 3 3]
% trim = for symmetric layouts, trim pd to only permutations
%        that result in unique ANOSIM or MDS (default = 1)
% max  = max size allowed for pd (default = 5000);
%        (if exceeded, NaN is returned)
%
% nPerms = total # of unique permutations of labels
% pd     = complete permutation distribution
%          (= NaN when nPerms > max)
%
% SEE ALSO: f_anosim

% -----Notes:-----
% This function is primarily for generating complete permutation
% distributions of grouping labels for symmetric distance matrices
% used as input for ANOSIM (f_anosim).

% -----Authors:---
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
%
% Subfunction multicomb written by special request
% by John D'Errico<derrico@flare.net>, 27-Mar-2002
%
% Thanks to Bob Clarke for explaining creation of unique
% permutation distributions when group size is symmetric

if (nargin < 2), trim = 1;    end; % trim permutation by default
if (nargin < 3), max  = 5000; end; % default max size = 5000

%-----Get size of complete permutation distribution:-----
grpsVar = unique(grps);    % extract unique grouping variables
noGrps  = length(grpsVar); % # of different groups
noLabels =length(grps);    % total number of grouping labels

for i=1:noGrps 
   grpSize(i) = sum(logical(grps==grpsVar(i))); % # of members in each group
   denom(i)   = factorial(grpSize(i));          % setup for denominator below
end;

% -----Get # of unique permutations:-----
% if group sizes are symmetric, need to divide by factorial(#groups)
% since only want permutations of labels resulting in UNIQUE
% values of R in ANOSIM or MDS configurations (R.K.Clarke, pers. comm.)

if (noGrps == sum(grpSize/grpSize(i)))
   symmVar = noGrps; % symmetric layout
else
   symmVar = 1;      % asymmetric layout
end;

% Get size of complete permutation distribution:
nPerms = factorial(noLabels)/(prod(denom)*factorial(symmVar));

% Generate permutation distribution if < max:
if nPerms <= max
   [p,pind] = multicomb(grpSize);
   pd = sortrows(pind);      % only want to keep the indices
   if (symmVar>1) & (trim>0) % layout is symmetric & trim is true
      nTotal = size(pd,1);   % total size of untrimmed permutation
      pd = pd([1:(nTotal/nPerms):nTotal],:); % take only 'unique' rows
   end;
else
   pd = NaN;
end;

% Experimental checking of results:
if (isnan(pd)==0) & (size(pd,1) ~= nPerms)
   error('nPerms and pd are NOT compatible!');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUBFUNCTION MULTICOMB %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,pind] = multicomb(k)
% -----Arguments:-----
%   k    = vector of # of members of each group
%   p    = combinations array
%   pind = indices

% ensure k is a row vector
k=k(:)';

% handle simple cases
if (sum(k)==0)
   % no objects
   p=[];
   if nargout>1
      pind=[];
   end
   return
elseif (sum(k>0)==1)
   % only one type of object left
   i=find(k);
   p=i*ones(1,k(i));
   if nargout>1
      pind=1:k(i);
   end
   return
end

% there are at least two different groups left
i=find(k);
p=[];
for j=i
   ki=k;
   ki(j)=ki(j)-1;
   temp=multicomb(ki);
   nt=size(temp,1);
   p=[p;[j*ones(nt,1),temp]];
end

% now convert the combinations to indices
% only do it if requested though. This ensures
% they are only computed at the very end
if nargout>1
   nk=length(k);
   ncomb=size(p,1);
   sumk=sum(k);
   pind=zeros(ncomb,sumk);
   csk=cumsum(k);
   for i=1:nk
      [c,r]=find(p'==i);
      cind=repmat(((1:k(i))+csk(i)-k(i))',1,ncomb);
      ind=sub2ind([ncomb,sumk],r(:),cind(:));
      pind(ind)=c(:);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%