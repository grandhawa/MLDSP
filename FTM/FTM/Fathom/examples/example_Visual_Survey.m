% Example of estimating the abundance fish within length classes for visual
% survey data
% 
% by David L. Jones, Aug-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             File: '.../examples/visual_survey.mat'           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The data in file 'visual_survey.mat' consist of observations of the
% abundance of fish encountered during visual survey transects along
% mangrove shoreline habitats in Biscayne Bay, Florida. See Jones et al.,
% 2010 for details.
% 
% Jones, D. L., J. F. Walter, E. N. Brooks, and J. E. Serafy. 2010.
%  Connectivity through ontogeny: fish population linkages among mangrove
%  and coral reef habitats. Mar. Ecol. Prog. Ser. 401: 245-258.

% Load the data:
load visual_survey.mat;

% -----Variables:-----
% fish = structure of 'transects x species' with the following fields:
%  .tra   = visual survey transect number
%  .id    = cell array specifying taxonomic code
%  .num   = # fish per 60 m^2 (rows = transect, cols = species);
%  .aveTL = observed average total length (mm TL)
%  .minTL = observed minmum total length (mm TL)
%  .maxTL = observed maximum length (mm TL)

% These data consist of the abundance and mean/minimum/maximum lengths of 10
% species of fish observed within each of 981 visual survey transects.
% For the first species (i = 1), estimate the number of fish within each of 5
% length classes for each of the 981 transects:

i      = 1; % specify first species
NaL_01 = f_numAtLength(fish.num(:,i), fish.aveTL(:,i), fish.minTL(:,i),fish.maxTL(:,i),1,0)
   
% NaL_01 = 
% 
%     L: [981x5 double] <- length class (mm TL)
%     A: [981x5 double] <- corresponding abundance (# fish per 60 m^2)
