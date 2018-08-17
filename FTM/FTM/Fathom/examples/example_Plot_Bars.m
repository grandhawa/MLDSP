% Examples of creating bar plots
% 
% by David L. Jones, Feb-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Create some data:
x      = [1:100]';
grp    = [ones(25,1); ones(25,1)*2; ones(25,1)*3; ones(25,1)*4];
labels = {'one' 'two' 'three' 'four'};

% Create bar plot for 1 variable:
figure;
f_plotBarsH(x,grp,0,labels,2);


% Create bar plot for multiple variables:
x2 = [x flipud(x)];
figure;
f_plotBarsV(x2,grp,0,labels,2);
