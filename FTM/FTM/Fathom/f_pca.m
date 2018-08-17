function result = f_pca(X,diag,method,tol,vRot,LS,ctr)
% - Principal Component Analysis of a data matrix
%
% Usage: result = f_pca(X,diag,method,tol,vRot,LS,ctr)
%
% X      = data matrix (rows = objects, cols = variables)          (default = 0)
% diag   = make diagnostic plots                                   (default = 0)
% method = use covariance (= 1), correlation (= 2), cross-products (= 3),
%          or Euclidean similarity (= 4) association matrix        (default = 1)
% tol    = trim eigenvalues < tol                                 (default = [])
% vRot   = use VARIMAX rotation                                    (default = 0)
% LS     = output Least-Squares scores                             (default = 0)
% ctr    = center data on column mean                              (default = 1)
%
% result = structure of results with the following fields:
%  .scores = coordinates of objects on each PC axis
%  .evects = eigenvectors (= loadings)
%  .evals  = eigenvalues for each PC axis
%  .expl   = percent and cumulative variance explained
%
% SEE ALSO: f_pcaPlot, f_pcoa

% -----Notes:-----
% This function is used to perform PCA using either the covariance
% or the correlation matrix of the input data. Use the COVARIANCE MATRIX
% when the variables are of the same kind, type, and scale; otherwise,
% use the CORRELATION MATRIX.
%
% Principal Component SCORES are the coordinates of the OBJECTS from
% the input data matrix in the new space defined by the Principal Component
% Axes. The EIGENVECTORS are the loadings of the original VARIABLES which, when
% multiplied by the original data, yield the Principal Component SCORES. The
% matrix of EIGENVECTORS has a row for each VARIABLE of the original data matrix
% and a column for each Principal Component Axis. The EIGENVALUES give the
% variance of the SCORES along each Principal Component Axis.
%
% The SCREE PLOT is a graphical method of evaluating how may PC Axes you need
% to retain in order to adequately represent the variation in the original
% data matrix--it provides the same information as EXPL.
%
% For some applications, you may want VARIMAX rotation + Least-Squares
% scores + data that are NOT centered (see Elmore & Richman, 2001)

% -----Details:-----
% This function implements PCA via Eigenanalysis of an association matrix formed
% from the original data. This matrix is decomposed into object space (U),
% variable space (V), and a diagonal matrix with singular variables along its
% diagonal (D) such that [association matrix] = [U*D*V], where eigenvectors = U
% and eigenvalues = the diagonal elements along D.

% -----References:-----
% Elmore, K. L., and M. B. Richman. 2001. Euclidean distance as a
%   similarity metric for principal component analysis. Monthly Weather
%   Review 129: 540-549.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, March-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: now based on f_eig vs. svd; origin of plot is marked
% Mar-2008: changed num2Str to num2str; skip plotting when # axes < 1, added
%           'tol'; white backgrounds in figures
% Mar-2010: renamed pflag to plt for consistency
% Feb-2011: removed redundant variable Xctr; added support for
%           cross-products and Euclidean similarity association matrices,
%           and least-squares scores
% Jun-2012: fixed scaling and plot title related to LS
% Dec-2012: results returned as a structure; plotting routines moved to a
%           separate function (f_pcaPlot)
% Jan-2013: default no diagnostic plot

% -----Set defaults & check input:-----
if (nargin < 2), diag   = 0;  end % default no diagnostic plot
if (nargin < 3), method = 1;  end % use covariance association matrix by default
if (nargin < 4), tol    = []; end % default don't trim eigenvalues
if (nargin < 5), vRot   = 0;  end % default no VARIMAX rotation
if (nargin < 6), LS     = 0;  end % default no least-squares scores
if (nargin < 7), ctr    = 1;  end % default center data on column mean

[nRows,nCols] = size(X);

if (nCols == 1), X(nRows,2) = 0; end; % add col of 0's if a col vector

if (LS>0) && (vRot<1)
   error('Least-squares scores (LS = 1) requires a VARIMAX rotation (vRot = 1)!');
end
% -------------------------------------

% Center data on column mean:
if (ctr>0), X = f_center(X); end

% Create association matrix:
switch method
   case 1
      A = 1/(nRows-1)*(X'*X); % covariance matrix (eq. 4.6 page 392)
   case 2
      A = corrcoef(X);        % correlation matrix
   case 3
      A = (X'*X);             % cross-products (Elmore & Richman,2001)
   case 4
      A = f_dis(X','eucs')';
   otherwise
      error('Unknown method!')
end

% Eigenanalysis:
[evects,evals] = f_eig(A);
evals          = evals'; % row vector of eigenvalues (1 for each axis)

% Discard eigenvalues < tol:
if ~isempty(tol)
   idx    = find (abs(evals)>tol);
   evects = evects(:,idx);
   evals  = evals(idx);
end

% Optional VARIMAX rotation:
if (vRot>0)
   evects = rotatefactors(evects,'Method','varimax');
end

% Principal Component Scores:
if (LS<1)
   if (method==2)
      scores = f_stnd(X) * evects; % correlation     (eq. 9.13 page 407)
   else
      scores = X * evects;         % everything else (eq. 9.4 page 394)
   end
else % Least-squares scores:
   scores = X * evects * f_inv(evects'*evects); % (eq. 8 in Elmore & Richman, 2001)
end

% Variation Explained:
varExplained  = abs(evals)/sum(abs(evals))*100; % Percent Variance Explained
cvarExplained = cumsum(varExplained);           % Cumulative Percent Explained
expl          = [varExplained; cvarExplained];  % combine into 1 matrix

% Specify title for plotting function (f_pcaPlot):
switch method
   case 1
      if (LS<1)
         txt = 'PCA Scores using Covariance Matrix';
      else
         txt = 'Least-Square Scores using Covariance Matrix';
      end
   case 2
      if (LS<1)
         txt = 'PCA Scores using Correlation Matrix';
      else
         txt = 'Least-Squares Scores using Correlation Matrix';
      end
   case 3
      if (LS<1)
         txt = 'PCA Scores using Cross-Products Matrix';
      else
         txt = 'Least-Squares Scores using Cross-Products Matrix';
      end
   case 4
      if (LS<1)
         txt = 'PCA Scores using Euclidean Similarity';
      else
         txt = 'Least-Squares Scores using Euclidean Similarity';
      end
end


% -----Create optional diagnostic plot:-----
if (diag>0)
   % Scree Plot:
   figure;               % open new window
   set(gcf,'color','w'); % set bg color to white
   plot(0:size(expl(2,:),2),[0 expl(2,:)],'-r.');
   title('Scree Plot');
   xlabel('PC Axes');
   ylabel('Cumulative % Variation Explained');
   set(gca,'XTick',[1:size(expl(2,:),2)]);
   set(gca,'YTick',[0:10:100]);
   axis([0 size(expl(2,:),2)+1 0 100]);
end
% -------------------------------------------

% Wrap results up into a structure:
result.scores = scores; % coordinates of objects on each PC axis
result.evects = evects; % eigenvectors (= loadings)
result.evals  = evals'; % eigenvalues for each PC axis
result.expl   = expl';  % percent and cumulative variance explained
result.txt    = txt;    % title for plotting function
