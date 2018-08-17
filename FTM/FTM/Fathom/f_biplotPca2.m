function f_biplotPca2(evects,scale,offset,sLabels)
% - create 2-d PCA distance biplot
%
% USAGE: f_biplotPca2(evects,scale,offset,sLabels)
%
% evects  = eigenvectors from f_pca
% scale   = scaling factor for descriptor vectors (default = 1)
% offset  = label offset (default = 0);
% sLabels = cell array of vector labels (if empty, autocreate)
%           e.g., sLabels = {'sal' 'tmp' 'elev'};
%
% SEE ALSO: f_pca, f_biplotEnv2, f_biplotSpecies, f_vectorfit

% -----Author(s):-----
% by David L. Jones, Mar-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp.

% -----Notes:-----
% This function is used to create a PCA distance biplot. Vectors are
% plotted for each variable comprising the original data matrix for which
% PCA was performed. The direction each vector points indicates the direction
% of increase of each variable (i.e., gradient). Length of vectors indicate
% the relative contribution each variable contributes to the formation of
% the reduced space plotted (here that space defined by the 1st 2 Principal
% Component axes).

% Mar-2008: changed & to &&

nRows = size(evects,1);

% ----- Check input and set default values: -----
if (nargin < 2), scale   = 1; end; % set default scaling to 1
if (nargin < 3), offset  = 0; end; % default is no label offset
if (nargin < 4), sLabels = num2cell([1:nRows]); end; % default vector labels

% If labels are not cell arrays, try forcing them:
if iscell(sLabels)<1, sLabels = num2cell(sLabels); end;

biplot = evects*scale; % scaled biplot

hold on;

% ----- plot the vectors: -----
for j = 1:nRows
   thisVector = [0,0;biplot(j,1),biplot(j,2)];
   plot(thisVector(:,1),thisVector(:,2),'-k');
   
   if (biplot(j,1)>=0 && biplot(j,2)>=0)
      deltaX = offset;
      deltaY = offset;
   elseif (biplot(j,1)<0 && biplot(j,2)>=0)     
      deltaX = -1*offset;
      deltaY =    offset;
   elseif (biplot(j,1)<0 && biplot(j,2)<0)   
      deltaX = -1*offset;
      deltaY = -1*offset;
   elseif (biplot(j,1)>=0 && biplot(j,2)<0)
      deltaX =    offset;
      deltaY = -1*offset;
   end;
   
   h = text(thisVector(2,1)+(deltaX*0.5),thisVector(2,2)+(deltaY*0.5),sLabels(j)); % label vectors
   set(h,'HorizontalAlignment','center');
end;

hold off;