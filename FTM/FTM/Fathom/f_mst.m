function  mst = f_mst(dis,dat)
% - Minimum Spanning Tree from a symmetric distance matrix
%
% USAGE: mst = f_mst(dis,dat)
%
% dis  = symmetric distance matrix
% dat  = optional Euclidean coordinates of objects corresponding to DIS
%        (necessary if resulting MST is to be plotted by f_plotNeigh)
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
% SEE ALSO: f_mst_mex f_plotNeigh, f_delaunay, f_dnn, f_gabriel, f_relNeigh, f_eigenMaps

% -----Notes:-----
% This function uses Kruskal's algorithm to calculate a Minimum Spanning Tree
% for objects based on pair-wise distances specified in a symmetric distance
% matrix. It is particularly useful for checking or interpreting the results of
% an ordination diagram by overlying the connections among objects. MST's are
% also useful for creating eigenvector maps and the .tDis is included for
% creating PCNM's.

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
% Jan-2014: updated documentation; overhauled for increased speed & efficiency:
%           replaced nested FOR loop and calls to 'ismember' with f_isPresent &
%           f_isAbsent functions
% May-2014: remove connections between elements where distance is 0

% -----Set defaults & check input:-----
if (nargin < 2), dat = NaN; end % default no coordinates

if (f_issymdis(dis) == 0)
   error('Input DIS must be a square symmetric distance matrix');
end
% ------------------------------------

n = size(dis,2); % size of distance matrix

% Get indices of lower tridiagonal elements:
[r,c] = find(triu(ones(n,n))==0);

% Unwrap distances, sort ascending:
[vec,key] = sort(f_unwrap(dis,0),'ascend');

% Get list of corresponding row/col indices:
list = [r(key) c(key)];

% Remove elements with 0 distances:
list(vec==0,:) = [];
vec(vec==0)    = [];

% Preallocate & initialize:
con       = [list(1,:);nan(n-2,2)]; % shortest connection
len       = [vec(1,:); nan(n-2,1)]; % corresponding branch length
list(1,:) = [];                     % remove from remaining pool
vec(1)    = [];

% Add to MST if only 1 of 2 objects are already present:
for i = 2:n-1 % only need n-1 branches to connect n objects
   % Get index to 1st row of LIST where col 1 is present & col 2 is absent from
   % CON (or vice versa):
   idx = find(f_isPresent(list(:,1),con(~isnan(con))) .* f_isAbsent(list(:,2),con(~isnan(con))) ...
      + f_isAbsent(list(:,1),con(~isnan(con))) .* f_isPresent(list(:,2),con(~isnan(con))),1);
   
   con(i,:)    = list(idx,:); % add this branch to mst
   len(i)      = vec(idx);    % add corresponding len
   list(idx,:) = [];          % remove connection from remaining pool
   vec(idx)    = [];          % remove length from remaining pool
end

% Rescale distance to 0-100 scale:
sca = len/max(len)*100;

% Create a symmetrical binary connection matrix:
B      = zeros(n,n);                       % initialize
idx    = sub2ind([n n],con(:,1),con(:,2)); % convert subscript to indices
B(idx) = 1;                                % indicate pairs that are connected
B      = B + B';                           % fill upper tridiagonal

% Create truncated distance matrix:
tDis          = dis;      % make a copy
t             = max(len); % truncation distance
tDis(dis > t) = t*4;      % truncated distance matrix

% -----Wrap results up into a structure:-----
mst.dat  = dat;
mst.dis  = dis;
mst.tDis = tDis;
mst.con  = con;
mst.len  = len;
mst.sca  = sca;
mst.B    = B;
