% Examples of creating ANOVA design matrices from categorical data:
% 
% by David L. Jones, Mar-2008
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/categorical.mat'

% -----References:-----
% Legendre, P. and M. J. Anderson. 1999. Distance-based redundancy analysis:
%   testing multispecies responses in multifactorial ecological experiments.
%   Ecological Monographs 69(1): 1-24. (See Appendix C)
% Anderson, M. J. 2003. XMATRIX: a FORTRAN computer program for calculating
%   design matrices for terms in ANOVA designs in a linear model. Department of
%   Statistics, University of Auckland, New Zealand.

% Load data file:
load categorical.mat

% -----Create binary dummy codes:-----
pos_b   = f_dummy(pos);
shade_b = f_dummy(shade);
% 
% Display results (use NaN's as a separator):
[pos_b repmat(NaN,size(pos_b,1),1) shade_b]
% ans =
% 
%      1   NaN     1     0
%      1   NaN     1     0
%      1   NaN     1     0
%      1   NaN     1     0
%      1   NaN     0     1
%      1   NaN     0     1
%      1   NaN     0     1
%      1   NaN     0     1
%      1   NaN     0     0
%      1   NaN     0     0
%      1   NaN     0     0
%      1   NaN     0     0
%      0   NaN     1     0
%      0   NaN     1     0
%      0   NaN     1     0
%      0   NaN     1     0
%      0   NaN     0     1
%      0   NaN     0     1
%      0   NaN     0     1
%      0   NaN     0     1
%      0   NaN     0     0
%      0   NaN     0     0
%      0   NaN     0     0
%      0   NaN     0     0


% -----Create contrast codes:-----
pos_c   = f_xMatrix(pos);
shade_c = f_xMatrix(shade);
% 
% Display results (use NaN's as a separator):
[pos_c repmat(NaN,size(pos_c,1),1) shade_c]
% ans =
%      1   NaN     1     0
%      1   NaN     1     0
%      1   NaN     1     0
%      1   NaN     1     0
%      1   NaN    -1     1
%      1   NaN    -1     1
%      1   NaN    -1     1
%      1   NaN    -1     1
%      1   NaN     0    -1
%      1   NaN     0    -1
%      1   NaN     0    -1
%      1   NaN     0    -1
%     -1   NaN     1     0
%     -1   NaN     1     0
%     -1   NaN     1     0
%     -1   NaN     1     0
%     -1   NaN    -1     1
%     -1   NaN    -1     1
%     -1   NaN    -1     1
%     -1   NaN    -1     1
%     -1   NaN     0    -1
%     -1   NaN     0    -1
%     -1   NaN     0    -1
%     -1   NaN     0    -1


% -----Create orthogonal constrast codes:-----
pos_o   = f_xMatrix(pos,1);
shade_o = f_xMatrix(shade,1);
% 
% Display results (use NaN's as a separator):
[pos_o repmat(NaN,size(pos_o,1),1) shade_o]
% ans =
% 
%             1          NaN     -0.70711     -0.70711
%             1          NaN     -0.70711     -0.70711
%             1          NaN     -0.70711     -0.70711
%             1          NaN     -0.70711     -0.70711
%             1          NaN       1.4142            0
%             1          NaN       1.4142            0
%             1          NaN       1.4142            0
%             1          NaN       1.4142            0
%             1          NaN     -0.70711      0.70711
%             1          NaN     -0.70711      0.70711
%             1          NaN     -0.70711      0.70711
%             1          NaN     -0.70711      0.70711
%            -1          NaN     -0.70711     -0.70711
%            -1          NaN     -0.70711     -0.70711
%            -1          NaN     -0.70711     -0.70711
%            -1          NaN     -0.70711     -0.70711
%            -1          NaN       1.4142            0
%            -1          NaN       1.4142            0
%            -1          NaN       1.4142            0
%            -1          NaN       1.4142            0
%            -1          NaN     -0.70711      0.70711
%            -1          NaN     -0.70711      0.70711
%            -1          NaN     -0.70711      0.70711
%            -1          NaN     -0.70711      0.70711


% -----Create orthogonal Helmert contrasts for balanced data:-----
fac1_h = f_helmert(fac1);
fac2_h = f_helmert(fac2);
% 
% Display results (use NaN's as a separator):
[fac1_h repmat(NaN,size(fac1_h,1),1) fac2_h]
% 
% ans =
%      2     0   NaN     1
%      2     0   NaN     1
%      2     0   NaN     1
%      2     0   NaN    -1
%      2     0   NaN    -1
%      2     0   NaN    -1
%     -1     1   NaN     1
%     -1     1   NaN     1
%     -1     1   NaN     1
%     -1     1   NaN    -1
%     -1     1   NaN    -1
%     -1     1   NaN    -1
%     -1    -1   NaN     1
%     -1    -1   NaN     1
%     -1    -1   NaN     1
%     -1    -1   NaN    -1
%     -1    -1   NaN    -1
%     -1    -1   NaN    -1
% 
% Create an interaction term by multiplying the main factors:
fac12_h  = repmat(fac1_h,1,size(fac2_h,2)) .* repmat(fac2_h,1,size(fac1_h,2));
fac12_h
% 
%      2     0
%      2     0
%      2     0
%     -2     0
%     -2     0
%     -2     0
%     -1     1
%     -1     1
%     -1     1
%      1    -1
%      1    -1
%      1    -1
%     -1    -1
%     -1    -1
%     -1    -1
%      1     1
%      1     1
%      1     1


% -----Create a contrast-coded interaction term:-----
[null,int] = f_xMatrix([pos shade]);
int
% int =
% 
%      1     0
%      1     0
%      1     0
%      1     0
%     -1     1
%     -1     1
%     -1     1
%     -1     1
%      0    -1
%      0    -1
%      0    -1
%      0    -1
%     -1     0
%     -1     0
%     -1     0
%     -1     0
%      1    -1
%      1    -1
%      1    -1
%      1    -1
%      0     1
%      0     1
%      0     1
%      0     1


% -----Create a contrast-coded nested term:-----
AB = f_xMatrix([A B],0,1)
% AB =
% 
%      1     0     0     0     0     0     0     0     0
%      1     0     0     0     0     0     0     0     0
%     -1     1     0     0     0     0     0     0     0
%     -1     1     0     0     0     0     0     0     0
%      0    -1     1     0     0     0     0     0     0
%      0    -1     1     0     0     0     0     0     0
%      0     0    -1     0     0     0     0     0     0
%      0     0    -1     0     0     0     0     0     0
%      0     0     0     1     0     0     0     0     0
%      0     0     0     1     0     0     0     0     0
%      0     0     0    -1     1     0     0     0     0
%      0     0     0    -1     1     0     0     0     0
%      0     0     0     0    -1     1     0     0     0
%      0     0     0     0    -1     1     0     0     0
%      0     0     0     0     0    -1     0     0     0
%      0     0     0     0     0    -1     0     0     0
%      0     0     0     0     0     0     1     0     0
%      0     0     0     0     0     0     1     0     0
%      0     0     0     0     0     0    -1     1     0
%      0     0     0     0     0     0    -1     1     0
%      0     0     0     0     0     0     0    -1     1
%      0     0     0     0     0     0     0    -1     1
%      0     0     0     0     0     0     0     0    -1
%      0     0     0     0     0     0     0     0    -1

