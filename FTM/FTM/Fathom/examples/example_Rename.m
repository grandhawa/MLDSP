% Example of renaming variables across multiple files
% 
% by David L. Jones, Aug-2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% In this example, two files containing parsed LA-ICP-MS data have a
% series of variables corresponding to values obtained for NIST-612
% external standard reference material (SRM). However, when these data were
% originally created the variables were named 'nist_A' instead of
% 'nist612_A', for example. Since data in files created subsequently
% included both NIST-612 and NIST-614 SRM's, we need to recursively process
% these files and change, say, 'nist_A' to 'nist612_A', etc.

% Work on COPIES of the original files (so you can repeat this example):
copyfile('./misc/run_1/101209_original.mat','./misc/run_1/101209.mat');
copyfile('./misc/run_2/101216_original.mat','./misc/run_2/101216.mat');

% -> Note the dot-slash ('./') notation specifies the folder you're
%    currently logged into; the data files reside within sub-folders of
%    that parent folder. The 'f_renameBatch' function is flexible enough to
%    be used in a variety of ways, but the PRESENT example requires you to
%    be logged in to the parent folder before calling the 'f_renameBatch'
%    function, NOT any of its child folders (i.e., NOT '/misc/',
%    /misc/run_1/' or  /misc/run_2/').

% Create a cell array specifying the target files:
cFile(1) = {'./misc/run_1/101209.mat'};
cFile(2) = {'./misc/run_2/101216.mat'};

% Specify the regex match-and-replace patterns to change the names of the
% appropriate variables in the target files:
ptn = 'nist_';    % match pattern
rep = 'nist612_'; % replacement

% Make sure the current working directly is the appropriate parent folder:
pwd
% 
% ans =
% 
% /Users/djones/work/mWork/Fathom/examples

% Rename the variables according to the regex patterns:
f_renameBatch(cFile,ptn,rep,1,1)
% 
% Processing file: ./misc/run_1/101209.mat...
% Processing file: ./misc/run_2/101216.mat...
% ->...DONE!
