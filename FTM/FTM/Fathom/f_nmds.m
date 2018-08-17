function mds = f_nmds(dis,dim,initial,diag,iter,tol)
% - Nonmetric Multidimensional Scaling (nMDS)
%
% Usage: mds = f_nmds(dis,dim,initial,diag,iter,tol);
%
% -----Input:-----
% dis     = symmetric dissimilarity matrix
% dims    = number of dimensions of solution                       (default = 2)
% initial = initial config: 0 = random, 1 = PCoA, or user-supplied coordinates
%                                                                  (default = 0)
% diag    = make diagnostic plots                                  (default = 0)
% iter    = max # iterations                                     (default = 200)
% tol     = convergence tolerance of stress function          (default = 0.0001)
%
% mds = structure of results with the following fields:
%  .scores = coordinates of solution
%  .stress = final stress of solution
%  .dim    = # of dimensions of solution
%  .rsq    = rank Mantel statistic (fitted distance vs. original dissimilarities)
%
% SEE ALSO: f_nmdsPlot, f_pcoa, f_pca

% -----Notes:-----
% This is a rewrite of f_nmds that now calls TMW's mdscale routine instead of
% a complied C-code MEX file

% -----Author:-----
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: check that MDSCALE is available; adjust plot bounds; plot
%           jittered axis if dim=1
% Mar-2013: updated to work with edits to f_pca; plotting routines moved to a
%           separate function (f_nmdsPlot)
% Jul-2014: removed input variable 'labels'; updated documentation; removed
%            unused variable 'n'

%-----Check input and set default values:-----
if (nargin < 2), dim     = 2;      end % default 2 dimensions
if (nargin < 3), initial = 0;      end % default initial random configuration
if (nargin < 4), diag    = 0;      end % efault don't create diagnostic plots
if (nargin < 5), iter    = 200;    end % default max iterations of 200
if (nargin < 6), tol     = 0.0001; end % default convergence criterion
if (nargin < 7), labels  = num2cell([1:size(dis,1)]); end; % default labels

% Make sure mdscale is available:
if ~(exist('mdscale','file')==2)
   error('The function requires MDSCALE from TMW''s Statistics Toolbox!');
end

if (f_issymdis(dis) == 0)
   error('Input DIST must be square symmetric distance matrix');
end

% Set starting configuration:
if isequal(initial,0)     % RANDOM
   startVar     = 'random';
   replicatesVar = 10; % # times to restart with new initial configuration
   
elseif isequal(initial,1) % PCoA
   startVar      = 'cmdscale';
   replicatesVar = 1;
   
else                      % USER-SUPPLIED
   if size(dis,1)~=size(initial,1)
      error('DIS & INITIAL should have the same # or rows!');
   end
   startVar      = initial;
   replicatesVar = 1;
   dim           = []; % dimensionality is inferred from starting configuration
end

criterionVar  = 'stress'; % stress1

%---------------------------------------------

% Set extra options:
opt = statset('Display','final','MaxIter',iter,'TolFun',tol,'TolX',0.0001);

% Annotate the display output:
if isequal(initial,0)
   fprintf('\n---------------------------------------------------------\n')
   fprintf('Selecting min STRESS from %d random start configurations:\n',replicatesVar);
   fprintf('---------------------------------------------------------\n')
end

% Run the NMDS routine:
[scores,stress,D] = mdscale(dis,dim,'Criterion',criterionVar,'Start',startVar,...
   'Replicates',replicatesVar,'Options',opt);

% Annotate the display output:
if isequal(initial,0)
   fprintf('\nBest configuration has STRESS = %1.7f\n',stress);
   fprintf('---------------------------------------------------------\n')
end

% Varimax rotation:
pca    = f_pca(scores,0);
scores = pca.scores; % extract scores from structure
clear pca;

% Make sure scores represent the proper dimension:
if dim==1
   scores = scores(:,1); % PCA adds a column of zeros when dim=1 ?
end

% Mantel Statistic comparing distances to dissimilarities:
dist      = f_dis(scores,'euc');      % fitted distances in Euclidean space
mantelRsq = (f_mantel(dis,dist,1))^2; % rank-correlation

% Wrap results up into a structure:
mds.scores = scores;
mds.stress = stress;
mds.dim    = size(scores,2);
mds.rsq    = (mantelRsq * 100);

% -----Create optional diagnostic plot:-----
if (diag>0)
   
   % Unwrap lower tridiagonal:
   dist_vec = f_unwrap(dist); % distance
   dis_vec  = f_unwrap(dis);  % dissimilarity
   D_vec    = f_unwrap(D);    % disparity
   
   % Sort according to disparity:
   [nul,idx] = sortrows([D_vec(:),dis_vec(:)]); % [Disparities Dissimilarities]
   
   figure;
   set(gcf,'color','w');    % set bg color to white
   plot(dis_vec,dist_vec,'bo',dis_vec(idx),D_vec(idx),'r.-');
   grid on;
   axis square;
   axis tight;
   xlabel( 'Dissimilarities' );
   ylabel( 'Distances/Disparities');
   title( sprintf( 'stress=%1.4f', stress) ,'FontSize',8);
   legend({'Distances' 'Disparities'},'Location','NW');
end
% ------------------------------------------
