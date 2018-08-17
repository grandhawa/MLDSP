function result = f_simper(X,grp,verb,xLabels)
% - similarity percentages (SIMPER) & species contributions among 2 groups
%
% USAGE: result = f_simper(X,grp,verb,xLabels);
%
% X     = biotic input matrix (rows = sites, cols = species abundances)
% grp   = column vector of integers specifying group memberhip of sites
%         e.g., grps = [1 1 2 2 1 1 2]';
% verb    = optionally send result to display  (default = 1)
% xLabels = cell array of species labels; if empty, autocreate
%            e.g., xLabels = {'sp1' 'sp2' 'sp3'};
% 
% result = structure with the following fields:
%  .avDisTot = total average among-group dissimilarity
%  .avDis    = average among-group dissimilarity by species
%  .avDis_SD = discrimination power of each species
%  .expl     = percent contribution to avDisTot
% 
% SEE ALSO: f_anosim, f_dis

% -----Notes:-----
% This function augments the results of an Analysis of Similarity (ANOSIM)
% performed using the f_anosim function. Since SIMPER is based on the
% Bray-Curtis coefficient, this function should be paired with f_anosim only
% when the latter was used with a Bray-Curtis dissimilarity matrix as input.
% Once you've determied there is a significant difference among groups using
% ANOSIM (global R), perform a series of pair-wise ANOSIM tests to identify
% which pairs of groups are significantly different. Then, use as input to this
% function a subset of your original (transformed) data corresponding to those
% observations making up 2 of the groups comprising ONE of the pairs found to be
% significantly different. Obviously, multiple calls to this function will be
% necessary if more than one of the pair-wise ANOSIM tests were deemed
% significant.
% 
% This function has been tested against Primer 5 for Windows and provides
% the same results.

% -----References:-----
% Clarke, K. R. and R. M. Warwick. 1994. Change in marine communities: an
% approach to statistical analysis and interpretation. Natural Environment
% Research Council, UK, 144 pp. [See page 7-3, eq. 7.1]

% -----Author:-----
% by David L. Jones, Apr-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


% -----Check input & set defaults:-----
if (nargin < 3), verb =  1; end % send output to display by default
if (nargin < 4), xLabels  = cellstr(num2str([1:size(X,2)]'))'; end % default X labels

grp = grp(:);    % force as a column vector
[n,p] = size(X); % get # rows, # variables
if size(grp,1) ~= n, error('X and GRP must have same # of rows!'); end

uGrp = unique(grp); % unique groups
nGrp = numel(uGrp); % # of unique groups
if nGrp ~= 2, error('GRP must code for only 2 groups (= 1 pair)!'); end

% If labels are not cell arrays, try forcing them:
if iscell(xLabels)<1, xLabels = num2cell(xLabels); end

% Make sure labels are of compatible size:
xLabels = xLabels(:); % force cell array into a column
if size(xLabels,1) ~= p
   error('Size of xLables doesn''t match # of X variables!')
end
% -------------------------------------
grp  = f_recode(grp); % recode groups to start at 1

idx1 = find(grp==1); % get index to rows of group 1
idx2 = find(grp==2); % get index to rows of group 2

nr_1 = numel(idx1);  % get # rows in group 1
blk  = numel(idx2);  % get # rows in group 2

D    = repmat(NaN,nr_1*blk,1); % preallocate
d    = repmat(NaN,nr_1*blk,p); % preallocate
idx  = [0 0];                  % initialize

% Modified version of Bray-Curtis, only inter- vs. intra-group site comparisons:
for i = 1:nr_1 % repeat for each row in group 1
   r        = repmat(X(idx1(i),:),blk,1);        % extract this row of group 1, replicate
   R        = X(idx2,:);                         % extract all rows of group 2
   idx      = idx(end)+1:idx(end)+blk;           % get index to this block of distances
   D(idx)   = sum(abs(r-R),2)./ sum(r+R,2);      % B-C dissimilarity among groups (all species)
   d(idx,:) = abs(r-R)./ repmat(sum(r+R,2),1,p); % dissimilarity separate by species
end

% Calculate statistics:
avDisTot  = mean(D);              % total average among-group dissimilarity
avDis     = mean(d);              % average among-group dissimilarity by species
SD        = std(d);               % standard deviation
avDis_SD  = avDis./SD;            % discrimination power of each species
expl      = repmat(NaN,p,1);      % preallocate
idx       = ~isnan(SD);           % prevent divide by NaN
expl(idx) = (avDis(idx)/avDisTot)*100; % percent contribution to avDisTot

% Wrap results up into a structure:
result.avDisTot = avDisTot; 
result.avDis    = avDis';
result.avDis_SD = avDis_SD';
result.expl     = expl;
 

% -----Send output to display:-----
if (verb>0)
   % Sort index for output:
   [expl_srt,idx] = sort(expl,1,'descend'); % sort by % contribution
   cum       = cumsum(expl_srt);            % get cumulative percentage
   
   fprintf('\n==================================================\n');
   fprintf('                    SIMPER:\n');
   fprintf('(Total Average Dissimilarity = %1.4f)\n',avDisTot);
   fprintf('--------------------------------------------------\n');
   j = 0; % initialize
   for i = idx'
      j = j+1;
      fprintf('%10.4f %10.4f %10.2f %%  %10.2f %%   %s\n',result.avDis(i),...
         result.avDis_SD(i),result.expl(i),cum(j),xLabels{i});
   end
   fprintf('--------------------------------------------------\n');
   fprintf('Col 1 = average among-group dissimilarity (avDis)\n');
   fprintf('Col 2 = discrimination power (avDis_SD)\n');
   fprintf('Col 3 = percent contribution (expl)\n');
   fprintf('Col 4 = cumulative percent contribution \n');
   fprintf('Col 5 = species\n');
end
% ---------------------------------
