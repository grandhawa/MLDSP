% Example of creating good uniform random numbers using Wichman & Hill's (2006)
% algorithm.
% 
% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Draw a sample of n=1000:
R = f_randWH(1000,1);
% f_pdf('randWH_1000'); % save as PDF

% Draw a larger sample to see how close it converges to uniform:
f_randWH(1000*100,1);
% f_pdf('randWH_large'); % save as PDF
% 
% -> looks like you need a really large sample size to get good uniform
% coverage.

% Compare this with Matlab's built-in random number generator:
figure; hist(rand(1000,1),25);
figure; hist(rand(1000*100,1),25);

% Provide your own random seed (for repeatability):
f_randWH(1000,1,[0.1 0.2 0.3 0.4]);

% Provide a randomly chosen seed:
rng('default'); % reset Matlab's random # generator (for repeatability)
seed = rand(4,1);
f_randWH(1000,1,seed);
