function rel = f_relNeigh(dat)
% - create a Relative Neighbor graph from 2-D spatial coordinates
%
% USAGE: rel = f_relNeigh(dat)
%
% dat = 2 column matrix specifying spatial coordinates
% 
% rel = structure of results with the following fields
%  .dat = original input coordinates
%  .dis = symmetric Euclidean distance matrix
%  .con = pairs of objects connected by triangle edges
%  .len = length of corresponding edges
%  .sca = edge length rescaled to 0-100
%  .B   = symmetric binary connection matrix; 1's indicate pairs of objects that
%         are connected, 0's indicate otherwise
%
% SEE ALSO: f_delaunay, f_dnn, f_gabriel, f_mst, f_eigenMaps

% -----Notes:-----
% Spatial coordinates should be in UTM (or similar) formats and not polar
% coordinates (i.e., Lat/Long).

% -----References:-----
% The main algorithm for creating Relative neighbors was ported to Matlab from
% Nicholas Lewin-Koh's C source code 'relative.C', which is used by the R
% function relneigh.R of the SPDEP package.

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (size(dat,2)>2)
   error('DAT must be a 2 column matrix of spatial coordinates!');
end
% -------------------------------------

% Extract x,y components:
x = dat(:,1);
y = dat(:,2);

% Preallocate/initialize:
n      = size(x,1);
con    = zeros(n*3,2);
no_rel = 1;           

% ----- Determine Relative Neighbors (after 'relative.C'):-----
for i = 1:n

   for j = (i+1):n
      rad = s_dist(x(i), y(i), x(j), y(j));

      for k = 1:n
         if ( (k~=i) && (k~=j) && (s_dist(x(i), y(i), x(k), y(k)) < rad) ...
               && (s_dist(x(j), y(j), x(k), y(k)) < rad) );
            break
         end
      end

      if (k == n)
         con(no_rel,1) = i;
         con(no_rel,2) = j;
         no_rel        = no_rel + 1;
      end

   end
end
% ---------------------------------------

% Remove empty connections:
con(~sum(con,2),:) = [];

% Make symmetrical:
con = unique([con;fliplr(con)],'rows');

% Create Euclidean distance matrix:
dis = f_euclid([x y]');

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
rel.dat = dat;
rel.dis = dis;
rel.con = con;
rel.len = len;
rel.sca = sca;
rel.B   = B;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 SUBROUTINES:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = s_dist(x1, y1, x2, y2)
% - calculate the distance between 2 points (after 'relative.C')
D = sqrt((x1-x2)^2 + (y1-y2)^2);
