function result = f_indVal(Y,grp,iter)
% - species indicator values
% 
% USAGE: result =  f_indVal(Y,grp,iter);
% 
% Y    = matrix of response variables (rows = obs, cols = species)
% grp  = column vector of integers specifying group membership for rows in Y
% iter = # iterations for permutation test                         (default = 0)
% 
% result = structure of results with the following fields:
%  .indVal   = species indicator power values (ranges from 0-100%) 
%  .p_indval = corresponding randomized probabilities
%  .IV       = all indicator values 
%  .p_IV     = corresponding randomized probabilities
%  .S        = Boolean values indicating groups each species displays its
%              maximum indicator value in (1=true, 0=false)
%  .G        = for each species (col), the actual group associated with each indVal

% -----Notes:-----
% A:      specificity is maximum when a species occurs in only 1 group
% B:      fidelity is maximum when a species occurs in all sites defining a group 
% indVal: the indicator value is maximum (= 100%) when a species occurs in
%         only one group and is present in all sites comprising that group.
% 
% An indicator species is a taxon that is characteristic of a group of
% sites because: (1) it mostly occurs only in one group; and (2) occurs in
% most of the sites comprising that group.
% 
% indVal = this metric is a measure of the magniture of a species ability to
% indicate (characterize) groups of sites, while the corresponding p-value
% is a measure of the statistical significance of that metric.
% 
% IV = this matrix simply returns all the indicator values used to
% determine maximum value used to determine indVal.
% 
% S = since indVal returns the maximum indicator value (IV) caluclated for each
% species, it would be useful to identify which group(s) yielded the
% maximum value for each species. The Boolean matrix S display this
% information in a matrix where rows correspond to groups and columns
% correspond to species. A value of 1 within a column identifies the group
% that yielded that maximum indicator value.
% 
% G = this matrix provides the actual group number associated with the
% maximum indicator value returned by indVal and the Boolean = 1 values
% returned in S. When species have duplicate maximum indicator values in more
% than one group, this matrix may be partially populated with NaN's that
% served as place holders. 
% 
% Note: Legendre & Legendre (2012) suggest interpreting p-values with
% caution if the same taxa are used to both: (1) define the groups and (2)
% caluclate species indicator values.

% -----References:-----
% Dufrene, M. and P. Legendre. 1997. Species assemblages and indicator
%   species: the need for a flexible asymmetrical approach. Ecol. Monogr.
%   67(3): 345-366.
% Legendre, P., and L. Legendre. 2012. Numerical ecology, 3rd English
%   edition. Developments in Environmental Modelling, Vol. 24. Elsevier
%   Science BV, Amsterdam. xiv + 990 pp. 

% -----Author:-----
% by David L. Jones, Apr-2013
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), iter = 0;  end % no permutation test by default

[nr,nc] = size(Y); % # obs, # variables

grp = grp(:); % force as column vector
if nr ~= size(grp,1), error('Y & GRP need same # of rows'); end
% -------------------------------------

% Process groups:
uGrp = f_unique(grp); % unique groups, unsorted
nGrp = length(uGrp);  % # unique groups
if (iter>0 && nGrp==1)
   error('There''s only 1 group, re-run with iter=0!')
end
   
% Calculate Specificity & Fidelity 
Yb = double(Y>0);   % convert to binary (presence/absence)
A  = nan(nGrp,nc); % preallocate
B  = A;
for i = 1:nGrp % repeat for each group
   idx    = find(grp==uGrp(i));                 % get index to rows of this group
   n      = numel(idx);                         % # sites in this group
   A(i,:) = mean(Y(idx,:),1);                   % mean abundance of this group
   B(i,:) = sum(Yb(idx,:),1) ./ repmat(n,1,nc); % relative frequency within group
end

% Convert to relative abundance across groups:
A = A ./ repmat(sum(A,1),nGrp,1);

% Calculate indicator values:
IV = (A .* B) * 100;

% Keep maximum values:
indVal = max(IV,[],1);

% Identify which group a species obtains its maximum indicator value from:
S = (IV - repmat(indVal,size(IV,1),1)) == 0;

% -----Show groups that correspond to species' indVal:-----
% Allow for multiple groups having the same maximum indicator value:
G = nan(size(S)); % initialize

for i = 1:size(S,2) % repeat for each species
   % Get index to group yielding the maximum IV:
   idx      = find(S(:,i) == 1);
   n        = numel(idx); % # groups for this species
   G(1:n,i) = idx;
end

% Trim rows of G that are all NaN's:
idx      = find( sum(isnan(G),2) > 0);
G(idx,:) = [];
% ---------------------------------------------------------


%-----Permutation tests:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   randIV(nGrp,nc,iter-1) = 0; % preallocate 3 dimensional results array
     
   for i = 1:(iter-1) % observed value is considered a permutation
      perm          = f_indVal(Y,f_shuffle(grp),0);   % permute obs (rows)
      randIV(:,:,i) = perm.IV;                        % permuted indicator values
   end
   j = ( (randIV - repmat(IV,[1 1 iter-1]) ) >= 0); % get randomized stats >= to observed statistic
   p_IV     = (sum(j,3)+1)./(iter);                 % count values & convert to probability
   p_indVal = sum(p_IV .* S);                         % p-values only for indVal's
else
   p_IV     = NaN;
   p_indVal = NaN;
end
%-----------------------------


% Wrap results up into a structure:
result.indVal   = indVal;
result.p_indVal = p_indVal;
result.IV       = IV;
result.p_IV     = p_IV;
result.S        = S;
result.G        = G;

