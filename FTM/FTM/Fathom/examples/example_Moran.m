% Examples of calculating Moran's coefficient of spatial autocorrelation
%
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/moran.mat'

% Reproduce the example provided for the 'moran' and 'moran.mc' functions in the
% 'SPDEP for R' package:
 
% Load the data file:
load moran.mat

% Calculate Moran's I for the crime data:
[MC,p] = f_moran(crime.y,crime.W,1000)
% MC = 0.51095
% p  = 0.001
% 
% The example from 'R' produced the following:
% 
% > sim1 <- moran.mc(crime, col.W, 999)
% > sim1
% 
% 	Monte-Carlo simulation of Moran's I
% 
% data:  crime 
% weights: col.W  
% number of simulations + 1: 1000 
%  
% statistic = 0.511, observed rank = 1000, p-value = 0.001
% alternative hypothesis: greater 
