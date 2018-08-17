function tri = f_delaunay(dat)
% - create a Delaunay triangulation from 2-D spatial coordinates
% 
% USAGE: tri = f_delaunay(dat)
% 
% dat = 2 column matrix specifying spatial coordinates
% 
% tri = structure of results with the following fields
%  .dat = original input data
%  .dis = symmetric Euclidean distance matrix
%  .con = pairs of objects connected by triangle edges
%  .len = length of corresponding edges
%  .sca = edge length rescaled to 0-100
%  .B   = symmetric binary connection matrix; 1's indicate pairs of objects that
%         are connected, 0's indicate otherwise
% 
% SEE ALSO: f_plotNeigh, f_dnn, f_mst, f_gabriel, f_relNeigh, f_eigenMaps

% -----Author:-----
% by David L. Jones, March-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2012: replaced f_euclid with f_dis

% -----Set defaults & check input:-----
if (size(dat,2)>2)
   error('DAT must be a 2 column matrix of spatial coordinates!');
end
% -------------------------------------

% Create Euclidean distance matrix:
dis = f_Dis(dat,'euc');
n   = size(dis,1);

% Get triangle vertices:
con = delaunayn(dat);

% Rearrange as pairs of edges:
con = [ con(:,1:2); con(:,2:3); [con(:,3) con(:,1)] ];

% Make sure it's completely symmetrical:
con = unique([con;fliplr(con)],'rows');

% Get indices to corresponding elements in distance matrix
idx = sub2ind([n,n], con(:,1), con(:,2));

% Create a symmetric binary connection matrix:
B      = zeros(n,n); % initialize
B(idx) = 1;
 
% Get connections after removing redundant edges:
BB                           = B;        % make a copy
BB(logical(triu(ones(n,n)))) = 0;        % force 0's in upper tri-diag
idx                          = find(BB); % get connections in lower tri-diag
con                          = [];       % re-initialize
[con(:,1) con(:,2)]          = ind2sub([n n],idx);

% Get corresponding edge lengths:
len = dis(idx);

% Rescale lengths to 0-100 scale:
sca = len/max(len)*100;

% -----Wrap results up into a structure:-----
tri.dat = dat;
tri.dis = dis;
tri.con = con;
tri.len = len;
tri.sca = sca;
tri.B   = B;
