function [result,bin] = f_disprof_clust(Y,dis,link,mc,iter,alpha,tol)
% - dissimilarity profile analysis (DISPROF) of a cluster analysis dendrogram
%
% USAGE: [result,bin] = f_disprof_clust(Y,'dis',link,mc,iter,alpha,tol)
%
% Y     = matrix of response variables (rows = obs, cols = variables)
% dis   = dissimilarity measure to apply to Y
%         (e.g., dis = 'bc'; see help for f_dis)
% link  = type of cluster analysis to perform:                     (default = 1)
%         1: 'average'  (UPGMA)
%         2: 'centroid' (UPGMC)
%         3: 'median'   (WPGMC)
%         4: 'single'   (Shortest Distance)
%         5: 'ward'     (Inner Squared Distance)
%         6: 'weighted' (WPGMA)
% mc    = progressive adjustment of p-values                       (default = 0)
% iter  = # permutation iterations for DISPROF analysis         (default = 1000)
% alpha = significance level                                    (default = 0.05)
% tol   = tolerance for rejecting a p-value > alpha            (default = 0.005)
%
% result = structure of results with the following fields:
%   .grp    = vector specifying group membership for each row of Y
%   .inc    = record of incremental changes made to GRP
%   .Pi     = Pi statisitic used to assess homogeneity of composition
%   .p      = corresponding randomized p-values
%   .p_idx  = corresponding rows of Y used in the assessment
%   .adj    = cell array of strings indicating: {'p-values adjusted'} or
%             {'p-values not adjusted'}
%   .Z      = cluster analysis tree linkages
%   .colS   = column vector indicating splits that were associated with
%             significant p-values (=1) and those that were not (=0).
%   .method = method used for constructing cluster analysis tree
% 
% bin = symmetric binary connectivity matrix defining pair-wise relationships as:
%        0: objects are members of the same cluster
%        1: objects are members of different clusters
%
% SEE ALSO: f_dis, f_disprof, f_disprof_clust_bin, f_disprof_clustPlot

% -----Notes:-----
% This function follows Clarke et al.'s (2008) description of similarity profile
% analysis (SIMPROF), but employs an equivalent approach based on dissimilarity
% profile analysis (DISPROF).
%
% This method is a form of 'agglomerative, hierarchical' cluster analysis.

% -----References:-----
% Clarke, K. R., P. J. Somerfield, and R. N. Gorley. 2008.  Testing null
%   hypotheses in exploratory community analyses: similarity profiles and
%   biota-environmental linkage. J. Exp. Mar. Biol. Ecol. 366:56-69.
% Legendre, P. & L. Legendre. 2012. Numerical ecology. 3rd English ed.
%   Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, Jan-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Jan-2014: improved format of output to display; now handles duplicate rows in
%            Y;
% Apr-2014: cv(cv==0) are now replaced with eps vs. [] to handle multiple sets
%           of duplicates; now optionally returns binary connectivity matrix;
%           removed 'complete linkage' option as it might form duplicates
% 
% Nov-2014: updated documentation

% -----Check input and set defaults:-----
if (nargin < 3), link  = 1;     end % default use UPGMA
if (nargin < 4), mc    = 0;     end % default don't adjust p-values
if (nargin < 5), iter  = 1000;  end % default iterations for permutation test
if (nargin < 6), alpha = 0.05;  end % default alpha of significance tests
if (nargin < 7), tol   = 0.005; end % default tol (p=0.055 is not significant)

% Check input for LINK:
if (~isequal(dis,'euc') && ismember(link,[2 3 5]))
   error('DIS must be ''euc'' when link = 2, 3, or 5!')
end

% Set the LINKTYPE:
switch link
   case 1
      linkTxt  = 'average';
      linkType = 'UPGMA';
   case 2
      linkTxt  = 'centroid';
      linkType = 'UPGMC';
   case 3
      linkTxt  = 'median';
      linkType = 'WPGMC';
   case 4
      linkTxt  = 'single';
      linkType = 'Shortest Distance';
   case 5
      linkTxt  = 'ward';
      linkType = 'Inner Squared Distance';
   case 6
      linkTxt  = 'weighted';
      linkType = 'WPGMA';
   otherwise
      error('LINK must be value from 1 to 6!');
end
% ---------------------------------------

nr = size(Y,1); % get # rows of input data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CLUSTER ANALYSIS TREE:                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform agglomerative hierarchical cluster analysis:
Z = linkage(f_unwrap(f_dis(Y,dis))',linkTxt);

% Extract cut-off values, descending order:
cv = flipud(Z(:,3));

% Replace any distances = zero with eps (in case duplicate rows occur in Y):
cv(cv==0) = eps;

% -----Determine cluster membership at each split:-----
cluT = nan(nr,nr-1); % initialize
for i = 1:numel(cv)
   % Get cluster number at each split:
   T = cluster(Z,'cutoff',cv(i),'criterion','distance');
   % Recode so cluster numbers are consistent with each call to 'cluster':
   cluT(:,i) = f_dummy2cat((f_dis(T,'euc') == 0));
end

% Recode so cluster numbers start with 1:
oldC = f_unique(cluT(:)); % old cluster number
nC   = numel(oldC);       % get # clusters
newC = 1:numel(oldC);     % new cluster number
clu  = nan(nr,nr-1);      % initialize
for i=1:nC
   clu(cluT==oldC(i)) = newC(i);
end
if ~isequal(clu(:,1),ones(size(Y,1),1))
   error('Column 1 of CLU doesn''t specify a SINGLE group!')
end

% Create binary version of 'clu' so 0's indicate clusters that haven't
% changed since the previous split (col) and 1's indicate those that have; this
% is used to avoid assessing the same cluster more than once by DISPROF:
cluN = clu;          % make a copy
uGrp = unique(cluN); % get unique groups
nGrp = numel(uGrp);  % get # groups
for i=1:nGrp
   % Get index to cols that show no change from the previous column:
   idxCol = find(sum(abs(diff((cluN==uGrp(i))')'))==0);
   idxCol = idxCol+1;
   for j = idxCol
      cluN(cluN(:,j)==uGrp(i),j) = 0; % flag with zeros
   end
end
cluB = (cluN>0); % create a binary version

% Cleanup:
clear i j oldC newC cluT cluN;
% -----------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Test of HOMOGENEITY OF COMPOSITION via DISPROF:               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc        = size(clu,2);  % get # cols
grp       = ones(nr,1);   % initalize as a single group
grpB      = ones(nr,1);   % initalize as a single group
inc       = nan(nr,nc-1);
Pi        = nan(nr,1);
p         = nan(nr,1);
p_idx{nr} = [];
cnt       = 0;            % initialize p-value counter
cntG      = 1;            % initialize group counter
colS      = zeros(1,nc);  % initialize significant col indicator

fprintf('\n=======================================================\n');
fprintf('            Cluster Analysis via DISPROF:  \n')
fprintf('-------------------------------------------------------\n');
fprintf('Evaluating the presence of...\n');
for j=1:(nc-1)
   % Update grouping vector with clusters defined in the next column:
   idx       = (grp>0);       % get index to rows not in terminal groups
   grp(idx)  = clu(idx,j);    % update rows not marked 'terminal'
   grpB(idx) = cluB(idx,j+1); % update corresponding binary vector
   inc(:,j)  = grp;           % record incremental changes to grouping vector
   
   % Get sorted list of unique groups, but flagged with grpB:
   uGrp = unique(grp.*grpB);
   
   % Remove those marked as terminal groups (<0) or haven't changed since the
   % previous column (= 0):
   uGrp(uGrp<=0) = [];
   
   if ~isempty(uGrp)
      nGrp = numel(uGrp); % get # groups (there should only be 1)
      
      for i=1:nGrp % Process each group:
         idxG = find(grp==uGrp(i)); % get index to rows of this group
         
         % Check if min # of rows are present:
         if (numel(idxG)<2)
            % Mark rows as terminal group with a negative number:
            grp(idxG) = grp(idxG) * -1;
         else
            % ------------------------------------------------------------------
            % Check for the presence of multivariate structure in this group:
            cnt = cnt + 1; % increment counter
            fprintf('  %d groups: ',cntG+1);
            disprof    = f_disprof(Y(idxG,:),dis,iter,0,0,0); % check rows for structure
            Pi(cnt)    = disprof.Pi;                          % collect Pi statistics
            p(cnt)     = disprof.p;                           % collect p-values
            p_idx{cnt} = idxG;                                % collect list of rows tested
            
            % --Optionally adjust p-values for multiple comparisions:--
            if (mc>0)
               p(cnt) = p(cnt)*cnt; % progressive Bonferroni method (L&L, 2012: p.745)
               if p(cnt)>1, p(cnt) = 1; end; % set values > 1 to 1
            end
            fprintf('Pi = %6.4f, p = %6.4f\n',Pi(cnt),p(cnt));
            
            % If no significant structure exists, no further splitting of group:
            if ((p(cnt) - alpha)>tol)
               % Mark rows as terminal group with a negative number:
               grp(idxG) = grp(idxG) * -1;
            else
               cntG    = cntG + 1;    % increment group counter
               colS(j) = colS(j) + 1; % indicate the split here was significant
            end
            % ------------------------------------------------------------------
         end
      end
   end
end
% ----------------------------------------------------------------------
fprintf('\n- # groups identified = %d \n',cntG);
fprintf('- # permutation iterations = %d \n',iter);
fprintf('- alpha = %6.4f \n',alpha);
if (mc>0)
   fprintf('- p-values adjusted for multiple testing: YES \n')
else
   fprintf('- p-values adjusted for multiple testing: NO \n')
end
fprintf('- dissimilarity metric = %s \n',dis);
fprintf('-------------------------------------------------------\n\n');

% Format for output (remove flag indicating terminal groups):
grp = abs(grp);
inc = abs(inc);

% Recode grouping vector as consecutive integers:
grpOld  = grp;              % make a copy
incOld  = inc;
grp     = f_recode(grp);    % recode as consecutive integers
uGrpOld = f_unique(grpOld); % get unsorted list of original groups
uGrp    = f_unique(grp);    % get unsorted list of new groups
nGrp    = numel(uGrpOld);   % get # of groups
for i=1:nGrp
   inc(incOld==uGrpOld(i)) = uGrp(i);
end

% Remove NaN's:
idxNaN        = isnan(p);
p(idxNaN)     = [];
p_idx(idxNaN) = [];

% Wrap results up into a structure:
result.grp    = grp;
result.inc    = inc;
result.Pi     = Pi;
result.p      = p;
result.p_idx  = p_idx;
if (mc>0)
   result.adj = {'p-values adjusted'};
else
   result.adj = {'p-values not adjusted'};
end
result.Z      = Z; % list of cluster linkages needed to draw dendrograms
result.colS   = colS;
result.method = linkType;

% Optionally return binary connectivity matrix:
% 0 = objects in same cluster, 1 = objects in different clusters
if (nargout >1),
   bin = ~(f_dis(grp,'euc')==0);
end
