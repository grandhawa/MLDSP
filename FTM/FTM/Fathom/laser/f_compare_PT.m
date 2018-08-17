function result = f_compare_PT(tra_A,tra_B)
% - compare two otolith profile transects
% 
% USAGE: result = f_compare(tra_A,tra_B);
%
% tra_A, tra_B = otolith profile transect data created by f_cps2ppm_PT,
%                f_extract_PT, or f_interpolate_PT
% 
% result = structure of results with the following fields:
%   .r   = correlation of corresponding masses among profiles
%   .dev = associated deviations among profiles
% 
% SEE ALSO: f_interpolate_PT, f_extract_PT, f_cps2ppm_PT

% -----Notes:-----
% This function compares two otolith profile transects after first checking to
% see they are of compatible size (use f_interpolate_PT first if they are not)

% -----Author:-----
% by David L. Jones, Aug-2013
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% -----Set defaults and check input:-----
% Extract data:
A.pos = tra_A.pos;
A.ppm = tra_A.ppm;
B.pos = tra_B.pos;
B.ppm = tra_B.ppm;
nc    = size(B.ppm,2); % get # columns in ppm

% Cleanup:
clear tra_A tra_B;

% Check if sizes are compatible:
if (~isequal(numel(A.pos),numel(B.pos)))
   error('TRA_A and TRA_B are not the same length, see f_interpolate_PT!');
end
% ---------------------------------------

% Calculate correlation and deviations of corresponding masses among profiles:
r = nan(1,nc); % preallocate
d = nan(1,nc);
for i = 1:nc
   r(i) = f_corr(A.ppm(:,i),B.ppm(:,i),0,0); % correlation
   d(i) = sum(abs(A.ppm(:,i) - B.ppm(:,i))); % sum of absolute deviation
end

% Wrap results up into a structure:
result.r   = r;
result.d   = d;
