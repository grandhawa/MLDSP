% Example of data reduction for LA-ICP-MS data
% 
% by David L. Jones, Aug-2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Import raw data:
raw_B = f_importXL('101209_B.csv');

% Save:
fname = 'exampleDataReduction';
saver;

% Select signal/background regions:
f_cpsParse(raw_B); % -> nist612_F to nist612_GG, macs_F-G, n =  5 otoliths

% Setup Internal Standard (IS):
IS.txt = {'Ca'};
IS.ppm = 40*10000; % Otolith

% Stack multiple STD's:
cSTD = {nist612_F nist612_FF nist612_G nist612_GG};

% Process data:
% dwell = 10 ms : dwell time of  the quadrupole
% adj   = 0     : adjust all values < LOD to 0
% spike = 'r'   : remove spikes using the Rosner Test
% drift = 1     : correct for instrument drift using linear interp.
% tol   = 0.55  : R2 must be at least 0.55 to use linear interp.
% verb  = 1     : show results in display
f_cps2ppm(o20_A,cSTD,'nist612',IS,10,0,'r',1,0.55,1);
% 
% ==============================================================================
%     'Analyte'    'ppm'          'LOD'          'mMole/Mole'    '# spikes'    'Drift R2'
%     'Li7'        [  0.29579]    [ 0.028704]    [ 0.0042698]    [       0]    [ 0.91308]
%     'Na23'       [   2648.6]    [   5.2702]    [    11.543]    [       0]    [ 0.84286]
%     'Mg24'       [   18.849]    [ 0.062083]    [  0.077704]    [       0]    [ 0.97218]
%     'P31'        [   125.39]    [   4.5297]    [    0.4056]    [       0]    [       0]
%     'Ca43'       [   400000]    [   35.904]    [      1000]    [       1]    [       0]
%     'Sc45'       [  0.11601]    [   0.0908]    [0.00025855]    [       4]    [       0]
%     'V51'        [ 0.047926]    [ 0.036853]    [9.4264e-05]    [       0]    [       0]
%     'Cr53'       [        0]    [  0.27981]    [0.00048098]    [       3]    [       0]
%     'Mn55'       [ 0.099284]    [ 0.073516]    [0.00018107]    [       1]    [       0]
%     'Fe57'       [   219.95]    [  0.91298]    [   0.39462]    [       0]    [ 0.84816]
%     'Co59'       [  0.22105]    [ 0.014159]    [0.00037582]    [       0]    [ 0.96357]
%     'Ni60'       [   0.3063]    [ 0.089473]    [0.00052288]    [       4]    [       0]
%     'Cu63'       [  0.70591]    [  0.13112]    [  0.001113]    [       0]    [ 0.75547]
%     'Zn64'       [   1.9096]    [  0.10614]    [  0.002926]    [       0]    [ 0.59565]
%     'Cu65'       [  0.86246]    [  0.14063]    [ 0.0013599]    [       0]    [ 0.59897]
%     'Ge72'       [  0.15059]    [  0.12002]    [0.00020772]    [       0]    [       0]
%     'Rb85'       [ 0.082167]    [ 0.017693]    [9.6325e-05]    [       2]    [ 0.80606]
%     'Sr88'       [   1880.6]    [  0.12195]    [    2.1505]    [       1]    [ 0.86216]
%     'Y89'        [        0]    [0.0045579]    [3.7911e-06]    [       0]    [ 0.99101]
%     'Cd114'      [ 0.082234]    [ 0.042774]    [7.3297e-05]    [       1]    [       0]
%     'Sn118'      [  0.50128]    [ 0.041494]    [ 0.0004231]    [       0]    [       0]
%     'Ba137'      [   6.3626]    [ 0.018038]    [ 0.0046422]    [       0]    [ 0.66995]
%     'Au197'      [        0]    [0.0093391]    [ 2.638e-06]    [       4]    [       0]
%     'Pb208'      [  0.34798]    [ 0.042558]    [0.00016827]    [       0]    [       0]
%     'Th232'      [        0]    [0.0029745]    [6.9857e-07]    [       1]    [ 0.84554]
%     'U238'       [0.0025237]    [0.0022237]    [1.0623e-06]    [       0]    [       0]
% 
% 
% ------------------------------------------------------------------------------
% SRM              = nist612
% LOD adjustment   = zero
% Spike removal    = Rosner Test
% Drift correction = linear
% NaN              = analyte not present in SRM
% 
% Drift R2 is the correlation of the linearly interpolated drift correction line:
%          0 = indicates Nearest Neighbor interpolation was used instead
%          - = indicates NO interpolation was performed
% ------------------------------------------------------------------------------