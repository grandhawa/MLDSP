% Example of processing otolith profile data
% by David L. Jones, Sep-2012

% -----Notes:-----
% Sampling parameters (Core-to-Edge profiles):
% Pre-ablation       = none
% Transect scans     = 64 um circle, 86% energy, 10 Hz, 10 um/sec
% Background         = 90 sec
% External Standards = NIST 612, 614, 610
% 
% Core-to-edge profile transects were performed on otoliths

% Raw LA-ICP-MS data have been previoulsy inported using the 'f_importXL'
% and 'f_cpsParse' functions.

% Load data:
load 120416_demo.mat

% Setup Internal Standard (IS):
IS.txt = {'Ca'};
IS.ppm = 40*10000;    % Otolith

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Process Profile Transect data:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust if < LOD  : no
% spike removal    : none
% drift correction : linear (R>0.55)
% filter/smooth    : no
% retain raw cps   : yes
tra_SM_23A = f_cps2ppm_PT(SM_PT_23A, {nist612_A nist612_B},10,'nist612',IS,10,[2],0,1,0.55,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plot profiles for selected analytes:   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify target masses:
tar = {'Na23' 'Sr88' 'V51' 'Ni60' };

% Plot filtered CPS data on a log-scale:
f_plot_PT(tra_SM_23A, tar,[],1,1,1);
