% Example of batch processing a set of raw, parsed otolith microchemistry
% spot data.
% 
% by David L. Jones, Jun-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% This file demonstrates how to batch process a set of raw LA-ICP-MS
% transient signal spot scan data that have already been previously
% imported by the 'f_importXL' function and parsed by the 'f_cpsParse'
% function. Following this example, separate data files should be placed in
% separate folders and the 'f_cps2ppm_SPOT' should be called from the
% parent folder.

% Create a cell array specifying files with parsed signal data to process:
cFile(1)  = {'./2010-12-09/101209.mat'};
cFile(2)  = {'./2010-12-16/101216.mat'};
% 
% -> Note the dot-slash ('./') notation specifies the folder you're
%    currently logged into; the data files reside within sub-folders of
%    that parent folder.

% Specify filename and save:
fname = 'all_data';
saver;

% Make sure Matlab is logged into the parent folder:
pwd
% ans =
% 
% /Users/djones/work/mWork/Fathom/laser/examples/exampleBatchProcess

% Make sure the data are organized in separate subfolders:
ls
% 2010-12-09			2010-12-16			exampleBatchProcessSpot.m

% Setup Internal Standard (IS):
IS.txt = {'Ca'}; IS.ppm = 40*10000;

% Process spot samples (just replicates A-C):
% adj   = 0    (set values <LOD = 0)
% spike = 'g'  (Grubbs Test)
% drift =  1   (drift correction via linear interpolation)
% tol   = 0.55 (R2 >= 0.55)
% verb  = 1    (verbose output)
% match = 1    (matches replicates A-C)
% day   = 1    (get date of acquistion from file)
% 
[oto,all] = f_cps2ppm_SPOT(cFile,'nist612',IS,10,0,'g',1,0.55,1,1,1);
% 
% Processing file ./2010-12-09/101209.mat...
%      o1003_A
%      o1003_B
%      o1003_C
%      o1004_A
%      o1004_B
%      o1004_C
%      ... 
% 
% Processing file ./2010-12-16/101216.mat...
%      o1062_A
%      o1062_B
%      o1062_C
%      o1063_A
%      o1063_B
%      o1063_C
%      ...
% 
% -> Note the above example uses 'match=1' whereby data from replicate spot
% scans (A-C) are averaged in the output.
