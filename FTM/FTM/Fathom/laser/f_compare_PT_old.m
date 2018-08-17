function [result,tra_Ai,tra_Bi] = f_compare_PT_old(tra_A,tra_B)
% - compare two otolith profile transects (legacy code)
% 
% USAGE: [result,tra_Ai,tra_Bi = f_compare_old(tra_A,tra_B);
%
% tra_A = otolith profile transect data created by the f_cps2ppm_PT function
% tra_B = ditto
% 
% result = structure of results with the following fields:
%   .r   = correlation of corresponding masses among profiles
%   .dev = associated deviations among profiles
% tra_Ai = interpolated values of tra_A (= NaN when no interpolation used)
% tra_Bi = interpolated values of tra_B (= NaN when no interpolation used)
% 
% SEE ALSO: f_compare_PT, f_interpolate_PT, f_extract_PT, f_cps2ppm_PT

% -----Notes:-----
% This functions compares two otolith profile transects after first checking to
% see they are of compatible size, and interpolating the smaller profile to
% match the larger when necessary.

% -----Author:-----
% by David L. Jones, Aug-2012
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

% -----Interpolate values when transects differ in size:-----
nA  = numel(A.pos); % get # dependent values
nB  = numel(B.pos);
int = 0; % initialize flag

if (nA<nB) % interpolate A
   int = 1;
   Ai.pos = linspace(A.pos(1),A.pos(end),nB);
   Ai.ppm = nan(nB,nc); % preallocate
   for i = 1:nc
      Ai.ppm(:,i) = interp1(A.pos,A.ppm(:,i),Ai.pos);
   end
   A = Ai;
   
elseif (nA>nB) % interpolate A
   int = 2;
   Bi.pos = linspace(B.pos(1),B.pos(end),nA);
   Bi.ppm = nan(nA,nc); % preallocate
   for i = 1:nc
      Bi.ppm(:,i) = interp1(B.pos,B.ppm(:,i),Bi.pos);
   end
   B = Bi;
end
% -----------------------------------------------------------

% Calculate correlation and devaiations of corresponding masses among profiles:
r = nan(1,nc); % preallocate
d = nan(1,nc);
for i = 1:nc
   r(i) = f_corr(A.ppm(:,i),B.ppm(:,i),0,0); % correlation
   d(i) = sum(abs(A.ppm(:,i) - B.ppm(:,i))); % sum of absolute deviation
end

% Wrap results up into a structure:
result.r   = r;
result.d   = d;
result.int = int;
switch int
   case 0
      tra_Ai = NaN;
      tra_Bi = NaN;
   case 1
      tra_Ai = A;
      tra_Bi = NaN;
   case 2
      tra_Ai = NaN;
      tra_Bi = B;
end

