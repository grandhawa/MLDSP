function clz
% - same as clear + clc + close all
% 
% SEE ALSO: cl

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2012: converted from a script to a function so 'MakeContents.m' will
%           include it in 'Contents.m'

% Clear variables from workspace:
evalin('caller','clear');

% Clear console:
clc;

% Close all figures:
close all;
