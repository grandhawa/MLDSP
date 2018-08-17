% Example of processing otolith core-to-edge profiles generated from (1)
% surface line scans and (2) line of spot scans
% 
% by David L. Jones, Jun-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% updated June-2012

% This file demonstrates how to process the raw transient signal data
% generated from LA-ICP-MS assays of otolith core-to-edge profiles. For
% this example, the raw data has already been previously imported by the
% 'f_importXL' function and parsed by the 'f_cpsParse' function. Note there
% is a specific sequence of loading and saving files so the f_cps2ppm_PS
% function can batch process the spot scans. There are also strict
% conventions for naming variables.

% -----Notes:-----
% Sampling parameters (Core-to-Edge profiles):
% Transect scans = 83 um circle, 86% energy, 10 um/sec, ~ 2 mm long
% Adjacent spots = 64 um circle, 86%, 60 sec duration
% Background     = 60 sec
% 
% PT06_A         = a single line scan
% PS06_01-31     = a series of 31 spot scans  
 
% Note: Prior to transect scans, the sampling region was preablated with a
% raster 2x of 150 um squares at 730 um/sec. Adjacent spot scans (64 um)
% were run inside of previous transect scan (83 um).

% Load the parsed raw data:
load microchem.mat;

% Save with a new file name:
fname = 'red_snapper';
saver;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Process Profile Transect data:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust if < LOD  : no
% spike removal    : none
% drift correction : linear (R>0.55)
% filter/smooth    : no
% retain raw cps   : yes
tra_06_A = f_cps2ppm_PT(PT06_A,{nist_A nist_AA nist_B nist_BB},10,'nist612',IS,10,2,0,1,0.55,0,0);

% -----Plot profile transect:-----
% 
% confidence envelopes : none
% LOD                  : no
% log scale            : yes
% filter/smooth data   : yes
% raw CPS data         : no
f_plot_PT(tra_06_A,{'V51' 'Ni60' 'Ba137'},[],0,1,1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Process Series of Spot Data:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pack external standard into a cell array (NOTE: the naming convention for
% this variable is very strict):
cSTD_PS06 = {nist_B nist_BB nist_C nist_CC nist_D nist_DD nist_E nist_EE nist_F nist_FF};

% Save file so 'cSTD_PS06' is written to disk and can be loaded by the
% f_cps2ppm_PS function:
saver;

% Process all spots in file 'red_snapper.mat' with the 'PS06' prefix: 
% 
% adjust if < LOD        : no
% spike removal          : Rosner test
% drift correction       : linear (R>0.55)
% verbose output         : yes
% get date from filename : no
spo_06 = f_cps2ppm_PS('red_snapper.mat','PS06','nist612',10,2,'r',1,0.55,1,0);

% -----Plot profile line-of-spots:-----
% LOD       : no
% log scale : yes
f_plot_PS(spo_06,{'V51' 'Ni60' 'Ba137'},0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Surface Transect vs. Series of Spot Scans      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_plot_PT(tra_06_A,{'Ba137'},[],0,1,1,0); % surface transect
f_plot_PS(spo_06,{'Ba137'},0,1);          % spot scans
theAxis = axis;
axis([min(spo_06.spot) max(spo_06.spot) theAxis(3:4)]);
