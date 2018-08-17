% Example off using Graphviz to visualize a minimal spanning tree (MST)
% by David L. Jones, Jun-2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    NOTES:                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% These examples require 'Graphviz' to be properly installed on your system. It
% is a cross-platform application that is freely available from
% http://www.graphviz.org.
% 
% Depending on your platform, you may also need to: (1) edit the
% 'f_graphviz_neato' function to specify the correct location of 'neato' defined
% in the 'pname' variable; and (2) edit the section that calls the
% external 'Preview.app' application if you are not on OS X. Windows users could
% call Adobe Distiller instead and Linux users could call Ghostscript.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          EXE Estuary Nematodes:                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'exe_nematodes.mat' consists of abundance data for 140 species of
% nematodes from 19 sites within the Exe estuary in the UK from Warwick (1971)
% and has the following variables:
% 
% bio.dat   = average counts of 140 spp of nematodes from 19 sites
% bio.txt   = cell array of corresponding column (species) labels
% bio.sites = cell array of site labels
% 
% Warwick, R. M. 1971. Nematode associations in the Exe estuary. J. Mar. Biol.
%  Assoc. U.K. 51: 439-454.

% Clear workspace:
clz;

% Load the data::
load exe_nematodes.mat bio;

% 4th-root transform the data:
Y = f_normal(bio.dat,'4');

% Create dissimilarity matrix:
dis = f_dis(Y,'bc');

% Create DOT file:
f_graphviz_mst(dis,'nematode_neato_all','all');

% Create undirected graph:
f_graphviz_neato('nematode_neato_all.dot')

% Create undirected graph based only on distances defined by the MST:
f_graphviz_mst(dis,'nematode_neato_mst','mst');
f_graphviz_neato('nematode_neato_mst.dot')

