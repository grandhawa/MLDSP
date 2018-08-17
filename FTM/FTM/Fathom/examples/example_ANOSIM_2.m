% Example of Two-Way ANOSIM with no replication
% 
% by David L. Jones, 2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/anosim2.mat'
 
% The file 'anosim2.mat' contains data concerning nematode abundances from
% Warwick (1971). There is 1 response variable 'dis' representing a
% Bray-Curtis symmetric distance matrix from 4th root transformed species
% abundances, and 2 factors: 'site' and 'time'.
% 
% Warwick, R. M. 1971. Nematode associations in the Exe estuary. J. Mar.
% Biol. Ass. U.K. 51:439?454.

% Load the data:
load anosim2.mat

% Perform the analysis:
[r,p] = f_anosim2(dis,site,time,1,1000);
