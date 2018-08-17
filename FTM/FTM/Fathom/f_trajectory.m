function [distance,steps] = f_trajectory(crds)
% - get distance along a trajectory
%
% USAGE: [length,steps] = f_trajectory(crds);
%
% crds     = matrix of coordinates in p-dimensional space
%
% distance = distance each point is from beginning of trajectory
% steps    = distance each point is from previous

% -----Author:-----
% by David L. Jones, Jul-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

dist     = f_dis(crds,'euc');    % symmetric Euclidean distance matrix
steps    = [0 diag(dist,-1)']';  % get diagonal below main diagonal
distance = cumsum(steps);        % compute cumulative distance

