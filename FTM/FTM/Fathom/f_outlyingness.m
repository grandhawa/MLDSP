function out = f_outlyingness(X)
% - normalized measure of outlyingness
% 
% USAGE: out = f_outlyingness(X)
% 
% X   = input data (rows = observations, cols = variables)
% out = outlyingness measures

% -----Notes:-----
% For each observation, an outlying measure is computed as n / sum(squared
% distance), which is then normalized by subtracting the median and divided by
% the median absolute deviation.
% 
% Values > 10 suggest an observation may be an outlier (Breiman & Cutler, 2003)

% -----References:-----
% Breiman, L., and A. Cutler. 2003. Manual on setting up, using, and
%  understanding Random Forests v4.0. Technical Report. 
%  ftp://ftp.stat.berkeley.edu/pub/users/breiman/Using_random_forests_v4.0.pdf

% -----Author:-----
% by David L. Jones, May-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Compute outlier measures:
SS = sum(f_dis(X,'euc').^2,2); % sum-of-squares for every observation

% Median & MAD:
med_SS = median(SS); % median
mad_SS = f_mad(SS);  % median absolute deviation

% Normalized Outlierness:
out = (SS - med_SS) / mad_SS;

% Set negative values to 0:
out(out<0) = 0;
