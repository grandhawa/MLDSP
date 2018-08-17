% Examples of 1-way ANOSIM
%  
% by David L. Jones, 2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/anosim.mat'

% The file anosim.mat contains data from Gray et al. (1990) and ships with
% Clarke's PRIMER program. There is 1 response variable "species" representing
% abundances of 174 species (columns) from 39 collection sites (rows), and
% 1 factor, grps, a row vector specifying group membership (i.e. categories
% representing distance from the center of an oilfield).
% 
% Gray, J. S., K. R. Clarke, R. M. Warwick, & G. Hobbs. 1990. Detection of
% initial effects of pollution on marine benthos: an example from the
% Ekofisk and Eldfisk oilfields, North Sea. Mar. Ecol. Prog. Ser.
% 66: 285-299.

load anosim.mat;                   % load the data
species2 = f_transform(species,1); % square-root transform species abundances
dis      = f_dis(species2,'bc');   % Bray-Curtis symmetric distance matrix 

[r,p] = f_anosim(dis,grps);
% ==================================================
%          1-way ANOSIM RESULTS:
% --------------------------------------------------
% Sorted Groupings:
%  1  1  1  1  1  1  2  2  2  2  2  2  2  2  2  2  3  3  3  3  3  3  3  3  3  3
%  3  3  4  4  4  4  4  4  4  4  4  4  4
% 
% Global Test:
%   R = 0.5413  p = 0.0010 (1000 of 408316439225392431104 possible perms) 
% 
% Pair-Wise Tests:
%   1  2: R = 0.5539  p = 0.0010 (1000 of 8008 possible perms) 
%   1  3: R = 0.8200  p = 0.0010 (1000 of 18564 possible perms) 
%   1  4: R = 0.9277  p = 0.0010 (1000 of 12376 possible perms) 
%   2  3: R = 0.1596  p = 0.0080 (1000 of 646646 possible perms) 
%   2  4: R = 0.7635  p = 0.0010 (1000 of 352716 possible perms) 
%   3  4: R = 0.5582  p = 0.0010 (1000 of 1352078 possible perms) 
% 
% ==================================================