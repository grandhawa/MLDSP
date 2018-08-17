function result = f_extract_PT(X,p1,p2)
% - extract subset of otolith core-to-edge profile transect
%
% USAGE: result = f_extract_PT(X,p1,p2);
%
% X  = structure of 'otolith profile' surface transect data created by f_cps2ppm_PT
% p1 = starting position of subset to extract (um)
% p2 = ending position of subset to extract   (um)
%
% result = structure of results with the following fields:
%   .txt     = cell array of text for each analyte element
%   .iso     = corresponding isotope number
%   .txt_iso = combined element + isotope labels
%   .s       = subset of original time interval                             (sec)
%   .speed   = speed of laser across sample surface                        (um/s)
%   .pos     = subset of relative position along transect                    (um)
%   .ppm     = subset of concentration of unknown; filtered/smoothed        (ppm)
%   .LOD     = subset of limits of detection                                (ppm)
%   .ratio   = subset of molar ratios to internal standard           (mMole/Mole)
%   .SRM     = name of SRM used for external calibration
%   .adj     = type of adjustment applied to values below LOD     (zero,LOD,none)
%   .spike   = type of spike removal applied to the time series
%   .drift   = method of drift correction
%   .tol     = tolerance for linear interpolation
%   .R2      = correlation of drift correction line
%   .sm      = tag indicating whether data were filter/smoothed
%   .cps     = subset of original raw cps data (when present)
%   .p       = starting/ending points of subset (um)
%   .idx     = index to corresponding elements of input structure
%
% SEE ALSO: f_cps2ppm_PT

% -----Author:-----
% by David L. Jones, Aug-2013
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input and set defaults:-----
% Check starting and ending points:
if ~(p1<p2)
   error('P1 must be less than P2!');
end

% Make sure subset is within bounds of original data:
minPos = min(X.pos);
maxPos = max(X.pos);
if ((p1<minPos) || (p1>maxPos))
   error('P1 does not fall within the bounds of X.pos!')
end
if ((p2<minPos) || (p2>maxPos))
   error('P2 does not fall within the bounds of X.pos!')
end
% ---------------------------------------

% Replace starting/ending points with observed value:
p   = [p1 p2]; % combine as 1 variable
idx = [0 0]; % preallocate
for i = 1:2
   xDiff  = abs(X.pos-p(i)); % difference between actual position and selected value;
   idxMin = find(xDiff == min(xDiff));
   idx(i) = max(idxMin);     % for ties, select the larger index
   p(i)   = X.pos(idx(i));
end

% Set indices of subset to extract:
idxS = idx(1):idx(2);

% Extract subset:
X.s     = X.s(idxS);
X.pos   = X.pos(idxS);
X.ppm   = X.ppm(idxS,:);
X.LOD   = X.LOD(idxS,:);
X.ratio = X.ratio(idxS,:);
if (isfield(X,'cps'))
   X.cps = X.cps(idxS,:);
end

% Wrap results up into a structure:
result     = X;            % rename for output:
result.p   = [p(1) p(2)] ; % starting/ending points of subset (um)
result.idx = idx;          % index to corresponding elements of input structure
