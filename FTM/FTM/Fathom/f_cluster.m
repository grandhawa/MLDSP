function idx = f_cluster(yDis,plt,yLabels)
% - UPGMA-based cluster analysis of a symmetric dissimilarity matrix
%
% USAGE: idx = f_cluster(yDis,plt,yLabels);
%
% yDis    = square symmetric dissimilarity matrix
% plt     = plot dendrogram                                        (default = 0)
% yLabels = cell array of corresponding labels, if empty autocreate
%            e.g., yLabels = {'A' 'B' 'C'};
%
% idx     = index to rows of yDis, sorted to correspond to ordering in dendrogram
%
% SEE ALSO: f_dis

% -----Author:-----
% by David L. Jones, Apr-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2013: fixed error with default yLabels

% -----Set defaults & check input:-----
if (nargin < 2), plt     = 0; end % default no dendrogram
if (nargin < 3), yLabels = cellstr(num2str([1:size(yDis,2)]'))'; end % default Y labels

if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix');
end

% Make sure labels are of compatible size:
yLabels = yLabels(:); % force cell array into a column
if size(yLabels,1) ~= size(yDis,1)
   error('yLabels & yDis don''t have the same # of rows!')
end
% -------------------------------------

% Unwrap lower tridiagonal of distance matrix:
vec = f_unwrap(yDis);

% Perform hierarchical cluster analysis using UPGMA:
Z = linkage(vec','average');

if (plt>0)
   figure;
   % set(gcf,'color','w'); % set background color to white
   
   [H,T,idx] = dendrogram(Z,0,'LABELS',yLabels,'ORIENTATION','right');
   
   % Arrange indices so they correspond with dendrogram:
   idx = flipud(idx(:));
   
   % Customize plot:
   box on;
   title('\bfUPGMA-based Cluster Analysis');
   xlabel('Dissimilarity');
else
   idx = NaN;
end

