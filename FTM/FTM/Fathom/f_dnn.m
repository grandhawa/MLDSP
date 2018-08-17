function dnn = f_dnn(dis,bnd,dat)
% - distance-based Nearest Neighbor Graph from 2-D spatial coordinates
%
% USAGE: dnn = f_dnn(dis,bnd,{dat})
%
% dis  = symmetric distance matrix
% bnd  = col vector specifying upper/lower bounds (e.g., bnd = [0 0.25])
% dat = optional Euclidean coordinates of objects corresponding to DIS
% 
% dnn = structure of results with the following fields
%  .dat = original input coordinates
%  .dis = symmetric Euclidean distance matrix
%  .con = pairs of objects connected by triangle edges
%  .len = length of corresponding edges
%  .sca = edge length rescaled to 0-100
%  .B   = symmetric binary connection matrix; 1's indicate pairs of objects that
%         are connected, 0's indicate otherwise
%
% SEE ALSO: f_plotNeigh, f_delaunay, f_gabriel, f_mst, f_relNeigh, f_eigenMaps

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 2), error('No BND has been specified!'), end;
if (nargin < 3), dat = NaN; end % default no coordinates

if (f_issymdis(dis) == 0)
   error('Input DIS must be a square symmetric distance matrix');
end

bnd = bnd(:); % force as a col vector

if bnd(1)>=bnd(2)
   error('Lower boundary in BND must be less than upper boundary!');
end

if (size(bnd,1)~=2)
   error('BND must have exactly 2 elements!')
end
% ------------------------------------

n = size(dis,1);
B = ((dis >= bnd(1)) & (dis<=bnd(2)));
B(logical(eye(n,n))) = 0; % force 0's along diagonal:

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
dnn.dat = dat;
dnn.dis = dis;
dnn.con = con;
dnn.len = len;
dnn.sca = sca;
dnn.B   = B;
