function [result,SRM] = f_cps2ppm(UNK,STD,SRM,IS,dwell,adj,spike,drift,tol,verb)
% - calculate analyte concentration from LA-ICP-MS transient signal data
%
% USAGE: [result,SRM] = f_cps2ppm(UNK,STD,'SRM',IS,dwell,adj,spike,drift,tol,verb);
%
% UNK = structure of measured SAMPLE signal
% STD = structure of measured STANDARD signal
%
% ...the preceeding structures have the following fields:
%  .s     = time interval (sec)
%  .cps   = count rate (counts per sec)
%  .txt   = cell array of text indicating measured elements
%  .iso   = corresponding isotope number
%  .hh    = hour of acquisition
%  .mm    = min of acquisition
%  .gDate = date of acquisition
%  .bg    = structure of corresponding background data (with the same fields)
%
% SRM = name of variable in file SRM.mat to use for Standard Reference Material
%       e.g., SRM = 'nist612'
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
%         (= 1) or nearest neighbor method (= 2)                   (default = 0)
% tol   = min value of R2 required to use linear vs. nearest neighbor
%         interpolation                                         (default = 0.55)
% verb  = optionally send results to display                       (default = 1)
%
% result = structure of results with the following fields:
%  .txt     = cell array of text for each analyte element
%  .iso     = corresponding isotope number
%  .txt_iso = combined element + isotope labels
%  .ppm     = concentration of unknown                                     (ppm)
%  .LOD     = limits of detection                                          (ppm)
%  .ratio   = molar ratios to internal standard                     (mMole/Mole)
%  .SRM     = name of SRM used for external calibration
%  .adj     = type of adjustment applied to values below LOD     (zero,LOD,none)
%  .spike   = type of spike removal applied to the time series
%  .nSpike  = total # of spikes removed
%  .drift   = method of drift correction
%  .tol     = tolerance for linear interpolation
%  .R2      = correlation of drift correction line
%  .bg      = mean background counts of UNK                                (cps)
%  .gDate   = cell array of Gregorian date of acquisition
% 
% SRM = structure of Standard Reference Material values used; provided to
%       allow calculation of 'percent relative error', etc.
%
% SEE ALSO: f_importXL, f_cpsParse, f_grubbs, f_rosner, f_spike,
%           f_cps2ppm_SPOT

% -----Notes:-----
% This function is used to process the laser ablation transient signal data of a
% single sample parsed into separate signal/background components by the
% 'f_cpsParse' function. Both UNK and STD inputs must have been previously
% created by the 'f_cpsParse' function. SRM identifies the reference material
% represented by the STD. IS identifies the analyte that serves as the internal
% standard and its concentration in the UNK (e.g., for an otolith IS.txt =
% {'Ca'}; IS.ppm = 40*10000).
%
% ADJ is used to optionally adjust values that fall below the limits of
% detection. Values that fall below LOD can be set to 0 (adj=0) if you consider
% these values cannot be distinguished from 0. On the other hand, setting these
% values to LOD (adj=1) provides a "maximum likely concentration" (Heinrich et
% al., 2003).
%
% SPIKE is used to optionally remove outliers from the transient signal of each
% analyte in both the UNK and STD. The Grubbs Test (spike = 'g') does this by
% estimating a critical value using alpha = 0.05. The Rosner Test examines
% the residuals between the original times series and one that has been
% smoothed via a Butterworth filter. Spikes can also be removed by
% identifying values that fall beyond a specific number of standard deviations 
% from the overall mean (e.g., spike = 2.5). Spikes identified by the
% Grubbs Test are replaced with values from a running mean, while those
% identified by the Rosner Test are replaced with values from the filtered
% time series.
%
% DRIFT is used to optionally correct for mass-specific changes in the
% instrument's sensitivity with time. The linear interpolation method (drift=1)
% corrects for drift by interpolating values of 'm' along a correction line
% created by regressing the time of acquisition ('t') of multiple STD's vs. the
% calculated value of 'm' for each STD, given the time of acquisition of the
% UNK. Use of a linearly interpolated drift correction line assumes there's a
% strong linear relationship between 'm' and 't', which may not be the case.
% Therefore if the correlation (R2) of 't' vs. 'm' is not above a specific
% tolerance (TOL) for a particlular analyte, nearest neighbor interpolation is
% used instead for that analyte. This function outputs the value of R2 for the
% linearly interpolated drift correction line of each analyte, though returns a
% value of 0 when the tolerance (TOL) value of R2 was not reached and the
% nearest neighbor method was used instead. The nearest neighbor method
% (drift=2) uses the value of 'm' associated with the STD acquired closest in
% time with that of the UNK. Correction for drift using either method requires
% that data from more than one STD be input to the function, which can be
% achieved by packing multiple STD's into a cell array (e.g., cSTD =
% {STD_01 STD_02 STD_03}).

% -----Misc:-----
% Normalized sensitivity (S) provides a cps/ppm measure and is used to convert
% count rates (or LOD) from cps to ppm.
%
% This function requires the files 'SRM.mat' and 'periodic_table.mat' be in
% your Matlab search path.
%
% This function has been tested against the 'AMS' java program and 'Pepita for
% Windows' and provides similar results.

% -----References:-----
% Halter, W., T. Pettke, C. A. Henrich, and B. Rothen-Rutishauser. 2002. Major
%   to trace element analysis of melt inclusions by laser-ablation ICP-MS:
%   methods of quantification. Chemical Oceanography 183:63-86
% Henrich, C. A., T. Petke, W. E. Halter, M. Aigner-Torres, A. Audetat, D.
%   Gunther, B. Hattendorf, D. Bleiner, M. Guillong, and I. Horn. 2003.
%   Quantitative multi-element analysis of minerals, fluid and melt inclusions
%   by laser-ablation inductively-coupled-plasma mass-spectrometry. Geochimica
%   et Cosmochimica Acta 67(18): 3473-3496.
% Jackson, S. E. 2008. Chapter 11: Calibration strategies for elemental
%   analysis by LA-ICP-MS. Pages 169?188 in Sylvester, P., editor. Laser
%   Ablation-ICP-MS in the Earth Sciences. Current practices and outstanding
%   issues. Mineralogical Association of Canada, Short Course Series Volume
%   40, Vancoouver, B.C.
% Longerich, H. P., S. E. Jackson, and D. Gunther. 1996. Laser ablation
%   inductively coupled pasma mass spectrometric transient signal data
%   acquisition and analyte concentration calculation. J. Analyt. Atom.
%   Spectrom. 11:899-904
% Longerich, H. P., S. E. Jackson, and D. Gunther. 1997. Laser ablation
%   inductively coupled pasma mass spectrometric transient signal data
%   acquisition and analyte concentration calculation: Errata. J. Analyt. Atom.
%   Spectrom. 12:391

% -----Author:-----
% by David L. Jones, Aug-2010
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Sep-2010: negative background-corrected count rates set to 0; removed R^2; SRM
%           is output with result to calculate 'percent relative error'
% Oct-2010: UNK_bg and STD_bg data are now brought in WITHIN the UNK and STD
%           structures as a separate field (*.bg); added acquisition time
%           (hh/mm); instrument drift correction; output more detailed; updated
%           documentation
% Nov-2010: combined element + isotope labels, drift correction now works
%           if an UNK is sampled AFTER all STD's (via 'extrap' in interp1)
% Apr-2011: added support for Rosner Test
% May-2011: nSpike = NaN when spike=0; force 'UNK.R = mean(UNK.cps,1)' to
%           produce a row vector for f_cps2ppmTS; acquisition time now set to
%           average (vs. time of first measurement)
% Jun-2011: now outputs mean background counts of UNK, S, and gDate;
%           updated documentation & error messages; final check for
%           negative concentrations and/or wrong standards were specified
% Oct-2011: re-arranged checking of user's input for drift and # STD's; use
%           'any' when checking for negative ppm's; improved adjustment of
%           ppm's below LOD
% Feb-2012: at some point earlier, this function no longer returns ".S"

% -----Check input & set defaults:-----
if (nargin <  6), adj   =  2;    end % default no adjustment if PPM < LOD
if (nargin <  7), spike =  0;    end % default no spike removal
if (nargin <  8), drift =  0;    end % default no drift correction
if (nargin <  9), tol   =  0.55; end % default tolerance for linear interpolation of m
if (nargin < 10), verb  =  1;    end % send output to display by default

% Check tolerance:
if (drift==1) && ((tol<=0) || (tol>=1))
   error('TOL must be in the range 0 < tol < 1!')
end
if (drift~=1)
   tol = NaN;
end

% If STD is not a cell array of multiple STD's, place it in one:
if ~iscell(STD)
   cSTD{1} = STD;
else
   f_rename STD cSTD; % rename if it was a cell array
end

% Get # of STD cell arrays
nCell = numel(cSTD);

% Check user's input for drift correction:
if (drift==0) && (nCell>1)
   error('You supplied multiple STD''s, but specified NO drift correction!');
end

% Check requirements for drift correction:
if (drift>0) && (nCell==1) 
   error('More than 1 STD is required to correct for instrument drift!');
end

% Make sure UNK includes a DATE field:
if ~ismember('gDate',fieldnames(UNK))
   gDate = {'NaN'};
else
   gDate = UNK.gDate;
end

% Make sure inputs include a BG field:
if ~ismember('bg',fieldnames(UNK))
   error(['UNK does not contain a background region!'])
end
for i=1:nCell
   if ~ismember('bg',fieldnames(cSTD{i}))
      error(['STD ' num2str(i) ' does not contain a background region!'])
   end
end

% Extract background data:
UNK_bg         = UNK.bg;
cSTD_bg{nCell} = NaN; % preallocate
%
for i = 1:nCell
   cSTD_bg{i} = cSTD{i}.bg;
end

% Check size of input:
if size(IS.txt,2)~=1
   error('Only 1 Internal Standard is currently supported!');
end

% Make sure STD's and their BG's measure the same elements + isotopes:
for i=1:nCell
   if (~isequal(cSTD{i}.txt,cSTD_bg{i}.txt) || ~isequal(cSTD{i}.iso,cSTD_bg{i}.iso))
      error('STD and STD_bg do not include the same isotopes!')
   end
end

% Make sure multiple STD's measure the same elements + isotopes:
if (nCell>1)
   for i=2:nCell
      if (~isequal(cSTD{1}.txt,cSTD{i}.txt) || ~isequal(cSTD{1}.iso,cSTD{i}.iso))
         error('Multiple STD''s  do not include the same isotopes!')
      end
   end
end

% Make sure UNK & STD measure the same elments + isotopes:
if (~isequal(UNK.txt,UNK_bg.txt) || ~isequal(UNK.iso,UNK_bg.iso))
   error('UNK and UNK_bg do not include the same isotopes!')
end
for i=1:nCell
   if (~isequal(UNK.txt,cSTD{i}.txt) || ~isequal(UNK.iso,cSTD{i}.iso))
      error('UNK and STD do not include the same isotopes!')
   end
end

% Make sure signal/background regions don't overlap:
if ( (min(UNK_bg.s) >= min(UNK.s) && min(UNK_bg.s) <= max(UNK.s)) ||...
      (max(UNK_bg.s) >= min(UNK.s) && max(UNK_bg.s) <= max(UNK.s)) )
   error('Time interval of UNK & UNK_bg overlap!')
end
for i=1:nCell
   if ( (min(cSTD_bg{i}.s) >= min(cSTD{i}.s) && min(cSTD_bg{i}.s) <= max(cSTD{i}.s)) ||...
         (max(cSTD_bg{i}.s) >= min(cSTD{i}.s) && max(cSTD_bg{i}.s) <= max(cSTD{i}.s)) )
      error('Time interval of STD %d & its background overlap!',i);
   end
end

% Check that file 'SRM.mat' is in your path:
if (exist('SRM.mat','file')~=2)
   error('''SRM.mat'' is not in Matlab''s search path!')
end

% Check that the 'SRM' specified is present:
if (size(whos('-file','SRM.mat',SRM),1)<1)
   eval(['error(''The variable ' SRM ' is not present within SRM.mat!'')'])
end

% Load SRM:
SRMname = SRM;                            % copy variable name for output
eval(['load(''SRM.mat'', ''' SRM ''')']); % load specified variable
eval(['f_rename ' SRM ' SRM'])            % rename loaded variable to 'SRM'

% Trim/sort SRM to match UNK & STD
[null,LOC] = ismember(UNK.txt,SRM.txt);
SRM.txt    = SRM.txt(LOC(LOC>0));
SRM.ppm    = SRM.ppm(LOC(LOC>0));
clear null LOC;

% Make sure UNK & trimmed SRM measure the same elemnts:
if (~isequal(SRM.txt,UNK.txt))
   error('SRM and UNK do not include the same elements!')
end

% Make sure Internal Stardard is present in UNK (and SRM):
if (ismember(IS.txt,UNK.txt)==0)
   error('UNK does not include the Internal Standard (IS)!');
end

% Check that file 'periodic_table.mat' is in your path:
if (exist('periodic_table.mat','file')~=2)
   error('''periodic_table.mat'' is not in Matlab''s search path!')
end

% Load Periodic Table, trim/sort to match UNK & STD:
load periodic_table.mat elements; % must be in your Matlab path
[null,LOC]         = ismember(UNK.txt,elements.txt);
elements.txt       = elements.txt(LOC(LOC>0));
elements.name      = []; % not needed
elements.atomicNum = []; % not needed
elements.atomicWt  = elements.atomicWt(LOC(LOC>0)); % g/Mole
clear null LOC;

% Make sure UNK & trimmed Periodic Table match:
if (~isequal(elements.txt,UNK.txt))
   error('''periodic_table.mat'' and UNK do not include the same elements!')
end

% Set up dwell times:
dwell = dwell(:)';   % force row vector
if (numel(dwell)==1) % constant dwell
   dwell = repmat(dwell,1,size(UNK.cps,2));
else
   if (numel(dwell) ~= size(UNK.cps,2))
      error('Size mismatch between DWELL and UNK!')
   end
end

% Set acquisition time (= average time of measurement) to minutes after
% midnight:
UNK.t = UNK.hh*60 + UNK.mm + mean(UNK.s)/60;
for i=1:nCell
   cSTD{i}.t = cSTD{i}.hh*60 + cSTD{i}.mm + mean(cSTD{i}.s)/60;
end
% -------------------------------------

% Add combined element + isotope labels:
UNK.txt_iso = cellstr(strcat(UNK.txt',regexprep( cellstr(num2str(UNK.iso')),'\s', '') ))';

% Get indices, sizes:
idx     = find(ismember(UNK.txt,IS.txt)==1); % index of Internal Standard in UNK (& SRM)
[Na,nc] = size(UNK.cps);        % # sweeps in Unknown signal, # of analytes
Nb      = size(UNK_bg.cps,1);   % # sweeps in Unknown background
Ns      = repmat(NaN,nCell,1);  % preallocate
for i=1:nCell
   Ns(i) = size(cSTD{i}.cps,1); % # sweeps in Standard signal
end

% Background-corrected count rates:
UNK.cps            = UNK.cps - repmat(mean(UNK_bg.cps),Na,1); % Unknown
UNK.cps(UNK.cps<0) = 0; % set negative background-corrected count rates to 0
for i = 1:nCell
   cSTD{i}.cps = cSTD{i}.cps - repmat(mean(cSTD_bg{i}.cps),Ns(i),1); % Standard
   cSTD{i}.cps(cSTD{i}.cps<0) = 0; % set negative background-corrected count rates to 0
end

% -----Optional spike removal:-----
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

% Mean of background-corrected count rates:
UNK.R = mean(UNK.cps,1);          % Unknown  (row vector)
for i=1:nCell
   cSTD{i}.R = mean(cSTD{i}.cps); % Standard (row vector)
end

% Slope of the calibration curve (eq. 2 in Longerich et al. 1997:Errata; see
% also eq. 3 of Jackson, 2008):
m    = repmat(NaN,nCell,nc); % initialize
mInt = repmat(NaN,1,nc);
t    = repmat(NaN,nCell,1);

% Collect slope (m) and times (t) for each STD:
for i = 1:nCell
   m(i,:) = (cSTD{i}.R./SRM.ppm).*(repmat( SRM.ppm(idx)./cSTD{i}.R(idx) ,1,nc));
   t(i)   = cSTD{i}.t;
end

% Optional correction for Drift:
switch drift
   case 0     % No correction
      mInt     = m;
      R2       = repmat({'-'},1,nc);
      driftTxt = 'none';
   case {1,2} % Interpolate m:
      R2 = repmat(NaN,1,nc); % preallocate
      for i=1:nc % Check drift correction line for each analyte:
         if sum(isnan(m(:,i))) % Don't try to interpolate NaN values:
            mInt(i) = NaN;
            R2(i)   = NaN;
         else
            corrVar = f_corr(t,m(:,i))^2;
            if ((drift==1) && (corrVar>=tol))
               mInt(i) = interp1(t,m(:,i),UNK.t,'linear','extrap');
               R2(i)   = corrVar;
            else
               mInt(i) = interp1(t,m(:,i),UNK.t,'nearest','extrap');
               R2(i)   = 0;
            end
         end
      end
      if (drift==1)
         driftTxt = 'linear';
      else
         driftTxt = 'nearest';
      end
   otherwise
      error ('Drift correction only supports method 1 or 2!')
end

% Normalized sensitivity in 'cps per ppm' (Longerich et al., 1996):
% S = (STD.R./SRM.ppm).*(repmat((UNK.R(idx)/STD.R(idx))*(IS.ppm/SRM.ppm(idx)),1,nc));
S = repmat(UNK.R(idx)/IS.ppm,1,nc) .* mInt;

% Concentration of analyte for Unknown in 'ppm' (Longerich et al., 1996):
ppm = UNK.R ./ S;

% Molar ratio's (# of milliMoles of analyte per Mole of Internal Standard):
% 1 mg = 0.001 grams
% 1 kg = 1,000 grams
%
% ppm    = ug/g = mg/kg
% ppm/1000      = g/kg
% Atomic weight = g/Mole
%
% Analyte: ppm / (g/Mole)                     = mMol/Kg;
% IS: (ppm/1000) / (g/Mole) = g/kg / (g/Mole) = Mol/Kg;
% Analtye:IS molar ratio                      = mMol/Mol;
%
ratio = (ppm ./ elements.atomicWt) ./ ...
   repmat((ppm(idx)/1000)/elements.atomicWt(idx),1,nc);

% Limit of Detection in 'ppm' (Longerich et al., 1996):
%   -if standard deviation = 0, replace with a count rate equal to 1 ion
%    detected during a single sweep of the quadrupole (sensu SILLS program)
BGstdv      = std(UNK_bg.cps); % standard deviation of the background
idx         = find(BGstdv==0); % analytes with SD = 0
oneIon      = (1*1000)./dwell; % after Halter et al., 2002  eq. 10b
BGstdv(idx) = oneIon(idx);
LOD         = ((3*BGstdv) ./ S) .* repmat(sqrt(1/Nb + 1/Na),1,nc);

% Optionally adjust values below LOD:
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

% -----Perform final check on results:-----
% Check for negative ppm's:
if any(ppm(:)<0)
   fprintf('\n')
   fprintf('==============================================================================\n');
   fprintf('There are negative values of PPM...use results with caution!\n')
   if (nCell>1 && (UNK.t<min(t) || UNK.t>max(t)))
      fprintf('UNK was sampled outside the bounds of multiple STD...use results with caution!\n');
   end
   fprintf('==============================================================================\n');
end
% Check that correct Standards were used:
if (nCell>1 && (UNK.t<min(t) || UNK.t>max(t)))
   fprintf('\n')
   fprintf('==============================================================================\n');
   fprintf('UNK was sampled outside the bounds of multiple STD...use results with caution!\n');
   fprintf('==============================================================================\n');
end
% -----------------------------------------

% Wrap results up into a structure:
result.txt     = UNK.txt;          % cell array of text for each analyte element
result.iso     = UNK.iso;          % corresponding isotope number
result.txt_iso = UNK.txt_iso;      % combined element + isotope labels
result.ppm     = ppm;              % concentration of unknown                      (ppm)
result.LOD     = LOD;              % limits of detection                           (ppm)
result.ratio   = ratio;            % molar ratios to internal standard      (mMole/Mole)
result.SRM     = SRMname;          % name of SRM used for external calibration
result.adj     = adjTxt;           % type of adjustment applied to values below LOD
result.spike   = spikeTxt;         % type of spike removal applied to the time series
result.nSpike  = nSpike;           % total # of spikes removed
result.drift   = driftTxt;         % type of drift applied
result.tol     = tol;              % tolerance for linear interpolation
result.R2      = R2;               % correlation of instrument drift correction line
result.bg      = mean(UNK_bg.cps); % mean background counts of UNK                 (cps)
result.gDate   = gDate;            % date of acquisition

% -----Send output to display:-----
if (verb>0)
   table.txt_iso = result.txt_iso;
   table.ppm     = result.ppm;
   table.LOD     = result.LOD;
   table.ratio   = result.ratio;
   table.nSpike  = result.nSpike;
   table.R2      = result.R2;
   
   headerCell = {'Analyte','ppm','LOD','mMole/Mole','# spikes','Drift R2'}; % Set up column labels
   resultCell = f_struct2flat(table);       % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   
   fprintf('\n==============================================================================\n');
   disp(resultCell);
   fprintf('\n------------------------------------------------------------------------------\n');
   fprintf('SRM              = %s\n',SRMname);
   fprintf('LOD adjustment   = %s\n',adjTxt);
   fprintf('Spike removal    = %s\n',spikeTxt);
   fprintf('Drift correction = %s\n',driftTxt);
   fprintf('NaN              = analyte not present in SRM\n\n');
   disp('Drift R2 is the correlation of the linearly interpolated drift correction line:')
   disp('         0 = indicates Nearest Neighbor interpolation was used instead')
   disp('         - = indicates NO interpolation was performed')
   fprintf('------------------------------------------------------------------------------\n');
end
% ---------------------------------
