function  mst = f_mst_old(dis,dat)
% - legacy version of f_mst
%
% USAGE: mst = f_mst(dis,dat)
%
% dis  = symmetric distance matrix
% dat  = Euclidean coordinates of objects corresponding to DIS (necessary if
%        resulting MST is to be plotted by f_plotNeigh)
% 
% mst = structure of results with the following fields
%  .dat  = original input data
%  .dis  = symmetric Euclidean distance matrix
%  .tDis = truncated distance matrix
%  .con  = pairs of objects connected by triangle edges
%  .len  = length of corresponding edges
%  .sca  = edge length rescaled to 0-100
%  .B    = symmetric binary connection matrix; 1's indicate pairs of objects that
%          are connected, 0's indicate otherwise
% 
% SEE ALSO: f_plotNeigh, f_delaunay, f_dnn, f_gabriel, f_relNeigh, f_eigenMaps

% -----Notes:-----
% This function uses Kruskal's algorithm to calculate a Minimum Spanning Tree
% for objects based on pair-wise distances specified in a symmetric distance
% matrix. It is particularly useful for checking or interpreting the results of
% an ordination by overlying an MST on an ordination plot. MST's are also useful
% for creating eigenvector maps and the .tDis is included for creating PCNM's.

% -----Author:-----
% by David L. Jones, March-2003
%
% This file is part of the 'FATHOM Toolbox for Matlab' and
% is released under the GNU General Public License, version 2.

% Apr-2003: minor changes to plot options
% Mar-2008: updated docs; changed & to &&; set bg color; commented out fprintf
% Apr-2008: added B; renamed variables; output to a structure; removed plotting
%           (use f_plotNeigh); made compatible with other neighbor functions;
%           added .tDis for creating PCNM's
% Jan-2013: updated documentation

% -----Set defaults & check input:-----
if (nargin < 2), dat = NaN; end % default no coordinates

if (f_issymdis(dis) == 0)
   error('Input DIS must be a square symmetric distance matrix');
end
% ------------------------------------

n = size(dis,2); % size of distance matrix

% Get indices of lower tridiagonal elements:
[r,c] = find(triu(ones(n,n))==0);

% Unwrap distances with row/col indices, sort ascending:
list = sortrows([r c f_unwrap(dis,0)],3);

con = list(1,:); % initialize mst with shortest branch

for i = 2:n-1 % only need n-1 branches to connect n objects
   for j = 1:size(list,1) % cycle thru list of all branches
      % Add to MST if only 1 of 2 objects is already present:
      if (sum([ismember(list(j,1),con(:,1:2)) ismember(list(j,2),con(:,1:2))])==1);
         con(i,:)  = list(j,:); % add this branch to mst
         list(j,:) = [];        % remove branch from list
         break                  % go to next i
      end
   end
end


% Rescale distance to 0-100 scale:
con(:,4) = con(:,3)/max(con(:,3))*100;

% Extract output:
len = con(:,3);
sca = con(:,4);
con = con(:,1:2);

% Create a symmetrical binary connection matrix:
B      = zeros(n,n);                       % initialize
idx    = sub2ind([n n],con(:,1),con(:,2)); % convert subscript to indices
B(idx) = 1;                                % indicate pairs that are connected
B      = B + B';                           % fill upper tridiagonal

% Create truncated distance matrix:
tDis          = dis;         % make a copy
t             = max(len);    % truncation distance
tDis(dis > t) = t*4;         % truncated distance matrix

% -----Wrap results up into a structure:-----
mst.dat  = dat;
mst.dis  = dis;
mst.tDis = tDis;
mst.con  = con;
mst.len  = len;
mst.sca  = sca;
mst.B    = B;
