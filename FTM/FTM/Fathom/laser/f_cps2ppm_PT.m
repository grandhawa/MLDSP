function result = f_cps2ppm_PT(UNK,STD,speed,SRM,IS,dwell,adj,spike,drift,tol,sm,cps)
% - process 'otolith profile' surface transect data
%
% USAGE: result = f_cps2ppm_PT(UNK,STD,speed,'SRM',IS,dwell,adj,spike,drift,tol,sm,cps);
%
% UNK = structure of measured SAMPLE signal
% STD = structure of measured STANDARD signal
%
% ...the preceeding structures have the following fields:
%  .s   = time interval (sec)
%  .cps = count rate (counts per sec)
%  .txt = cell array of text indicating measured elements
%  .iso = corresponding isotope number
%  .hh  = hour of acquisition
%  .mm  = min of acquisition
%  .bg  = structure of corresponding background data (with the same fields)
%
% speed = translation speed of laser across sample surface                (um/s)
% SRM   = name of variable in file SRM.mat to use for Standard Reference Material
%         e.g., SRM = 'nist612'
%
% IS = structure of Internal Standard with the following fields:
%  .txt  = cell array of text indicating element
%  .ppm  = corresponding concentration in Unknown
%
% dwell = dwell time of the quadrupole (msec); use a scalar for constant
%         and a row vector for variable dwell times
% adj   = optionally adjust concentrations that fall below LOD to:
%         0 (adj=0), LOD (adj=1), or make no adjustment (adj=2, default)
% spike = optional spike removal of background-corrected signal data via Grubbs
%         Test (= 'g'), Rosner Test (='r'), or STDEV (= 1 to 4)    (default = 0)
% drift = optionally correct for instrument drift using linear interpolation
%         (= 1), nearest neighbor method (= 2) or none (= 0, default)
% tol   = min value of R2 required to use linear vs. nearest neighbor
%         interpolation                                         (default = 0.55)
% sm    = filter/smooth time-series data                           (default = 0)
% cps   = retain raw cps data                                      (default = 0)
%
% result = structure of results with the following fields:
%  .txt     = cell array of text for each analyte element
%  .iso     = corresponding isotope number
%  .txt_iso = combined element + isotope labels
%  .s       = original time interval                                       (sec)
%  .speed   = speed of laser across sample surface                        (um/s)
%  .pos     = relative position along transect                              (um)
%  .ppm     = concentration of unknown; filtered/smoothed                  (ppm)
%  .LOD     = limits of detection                                          (ppm)
%  .ratio   = molar ratios to internal standard                     (mMole/Mole)
%  .SRM     = name of SRM used for external calibration
%  .adj     = type of adjustment applied to values below LOD     (zero,LOD,none)
%  .spike   = type of spike removal applied to the time series
%  .drift   = method of drift correction
%  .tol     = tolerance for linear interpolation
%  .R2      = correlation of drift correction line
%  .sm      = tag indicating whether data were filter/smoothed
%  .cps     = original raw cps data
%
% SEE ALSO: f_plot_PT, f_export_PT, f_cps2ppm_PS, f_cps2ppm

% -----Notes:-----
% This function calculates analyte determinations for each observation of
% LA-ICP-MS transient signal data, then filters/smooths the data according
% to the method of Sinclair et al. 1998 (i.e., 11-point running median/average).
%
% Care must be taken when adjusting concentrations that fall below LOD's,
% because misleading results can arise if this is done before smoothing the
% data (and plotting smoothed ppm's alongside unsmoothed LOD's).

% -----References:-----
% Sinclair, D. J., L. P. J. Kinsley, & M. T. McCulloch. 1998. High
%  resolution analysis of trace elements in corals by laser ablation ICP-MS.
%  Geochim. Cosmochim. Acta 62:1889-1901.

% -----Author:-----
% by David L. Jones, May-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2011: set default TOL to 0.55; updated error msg; .s now added to
%           UNK_temp; renamed from *TS to *_PT
% Oct-2011: updated documentation; fixed setting of defaults; now filters
%           LOD;
% Dec-2011: data filtering now optional
% Apr-2012: optionally retain original raw cps data because someone may want to
%           examine analyte cps to total cps ratios.

% -----Check input & set defaults:-----
if (nargin <  7),  adj    = 2;    end % default no adjustment if PPM < LOD
if (nargin <  8),  spike  = 0;    end % default no spike removal
if (nargin <  9),  drift  = 0;    end % default no drift correction
if (nargin <  10), tol    = 0.55; end % default tolerance for linear interpolation of m
if (nargin <  11), sm     = 0;    end % default don't filter/smooth times-series data
if (nargin <  12), cps    = 0;    end % default don't retain original raw cps data

verb  = 0; % don't send output to display
% -------------------------------------

[nr,nc] = size(UNK.cps); % get size of input

% If STD is not a cell array of multiple STD's, place it in one:
if ~iscell(STD)
   cSTD{1} = STD;
else
   f_rename STD cSTD; % rename if it was a cell array
   if (drift==0)
      error('You supplied multiple STD''s, but specified NO drift correction!');
   end
end

% -----Optionally retain original raw CPS data:-----
if (cps>0)
   raw_cps = UNK.cps;
else
   raw_cps = NaN;
end

% -----Optional spike removal:-----
nCell = numel(cSTD);  % get # of STD cell arrays
if (spike~=0)
   nSpike = 0; % initialize
   if isequal(lower(spike),'g') % Grubbs Test:
      [UNK.cps,idxS] = f_grubbs(UNK.cps,1,0.05);
      nSpike         = nSpike + sum(idxS); % count # of spikes
      for i = 1:nCell
         cSTD{i}.cps = f_grubbs(cSTD{i}.cps,1,0.05);
      end
      spikeTxt = 'Grubbs Test';
   elseif isequal(lower(spike),'r') % Rosner Test:
      [UNK.cps,idxS] = f_rosner(UNK.cps,1,0.5,0.1,2);
      nSpike         = nSpike + sum(idxS); % count # of spikes
      for i = 1:nCell
         cSTD{i}.cps = f_rosner(cSTD{i}.cps,1,0.5,0.1,2);
      end
      spikeTxt = 'Rosner Test';
   elseif (spike>=1 && spike<=4)
      [UNK.cps,idxS] = f_spike(UNK.cps,1,spike);
      nSpike         = nSpike + sum(idxS); % count # of spikes
      for i = 1:nCell
         cSTD{i}.cps = f_spike(cSTD{i}.cps,1,spike);
      end
      spikeTxt = [num2str(spike) ' STDEV'];
   else
      error('Spike must be ''g'', ''r'', or range from 1-4!')
   end
else
   nSpike   = repmat({'-'},1,nc);
   spikeTxt = 'none';
end
% ---------------------------------

UNK_temp = UNK;        % make a copy
ppm      = nan(nr,nc); % preallocate
LOD      = nan(nr,nc);
ratio    = nan(nr,nc);
R2       = nan(1,nc);

fprintf('\nProcessing %d observations of %d analytes:\n',nr,nc);

% Process each observation separately:
for i = 1:nr
   
   % Show progress in display:
   if (mod(i,100)==0)
      fprintf('Processing observation %d...\n',i);
   end
   UNK_temp.cps = UNK.cps(i,:);
   UNK_temp.s   = UNK.s(i,:);
   % No spike removal or adjust ppm's < LOD's here:
   temp         = f_cps2ppm(UNK_temp,cSTD,SRM,IS,dwell,2,0,drift,tol,verb);
   ppm(i,:)     = temp.ppm;
   LOD(i,:)     = temp.LOD;
   ratio(i,:)   = temp.ratio;
   
   % Collect R2 only once:
   if (i==1)
      R2 = temp.R2;
   end
end
fprintf('Done!\n');

% Convert acquisition time to position
pos = (UNK.s-min(UNK.s))*speed;

% Optionally filter/smooth times-series data:
if (sm>0)
   ppm    = f_filterSinclair(ppm);
   LOD    = f_filterSinclair(LOD);
   ratio  = f_filterSinclair(ratio);
   sm_txt = {'data filtered'};
else
   sm_txt = {'data not filtered'};
end

% Optionally adjust values below LOD (after filtering):
if (adj==0) || (adj==1)
   idx = find(ppm < LOD);     % index to values below LOD
   switch adj
      case 0
         ppm(idx) = 0;        % set values to 0
         adjTxt   = 'zero';
      case 1
         ppm(idx) = LOD(idx); % set values to LOD
         adjTxt   = 'LOD';
   end
else
   adjTxt = 'none';
end

% Wrap results up into a structure:
result.txt     = temp.txt;     % cell array of text for each analyte element
result.iso     = temp.iso;     % corresponding isotope number
result.txt_iso = temp.txt_iso; % combined element + isotope labels
result.s       = UNK.s;        % original time interval (sec)
result.speed   = speed;        % speed of laser across sample surface         (um/s)
result.pos     = pos;          % relative position along transect               (um)
result.ppm     = ppm;          % concentration of unknown                      (ppm)
result.LOD     = LOD;          % limits of detection                           (ppm)
result.ratio   = ratio;        % molar ratios to internal standard      (mMole/Mole)
result.SRM     = temp.SRM;     % name of SRM used for external calibration
result.adj     = adjTxt;       % type of adjustment applied to values below LOD
result.spike   = spikeTxt;     % type of spike removal applied to the time series
result.nspike  = nSpike;       % total # of spikes removed
result.drift   = temp.drift;   % type of drift applied
result.tol     = temp.tol;     % tolerance for linear interpolation
result.R2      = R2;           % correlation of instrument drift correction line
result.sm      = sm_txt;       % tag indicating whether data were filtered/smoohted
result.cps     = raw_cps;      % optionally retain original cps data
