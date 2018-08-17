function id = f_cdaClass(x,y,z)
% - classify unknown observations for f_cda
%
% USAGE: id = f_cdaClass(x,y,z)
%
% x  = input data having known group membership
% y  = column vector of integers specifying group membership
% z  = input data to classify
%
% id = identification of observations in z according to group
% 
% SEE ALSO: f_cda, f_cdaCV

% -----Author:-----
% by David L. Jones, Dec-2003
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.
%
% This code is provided as is, with no guarantees.

% -----Set defaults and check input:-----
if (nargin < 4), verb  = 1; end; % send output to display by default
if size(x,2) ~= size(z,2)
   error('X and Z must have same # columns !'); 
end
% --------------------------------------

n = size(z,1); % # of rows (observations) of unknowns

% CDA on known data:
result_x   = f_cda(x,y,1,0,0,0);
Cvects     = result_x.Cvects;
centroids  = result_x.centroids;

% preallocate:
classifiedGrp = repmat(NaN,n,1);

% Center Z observations using means of X:   
z_ctr = z - repmat(mean(x),size(z,1),1);

% Project unknown obs in canonical space:
scores_z = z_ctr*Cvects;

% Find group centroid closest to unknown obs
% and classify accordingly:      
for i=1:n
   dist = [f_euclid([scores_z(i,:);centroids]')].^2; % squared distances
   centroidDist = dist(2:end,1);
   classifiedGrp(i) = find(centroidDist == min(centroidDist));
end

id = classifiedGrp;

