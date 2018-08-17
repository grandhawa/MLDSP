function result = f_cps2ppm_PS(fname,UNKname,SRM,dwell,adj,spike,drift,tol,verb,day)
% - process 'otolith profile' data collected as a series of spots
%
% USAGE: result = f_cps2ppm_PS('fname','UNKname','SRM',dwell,adj,spike,drift,tol,verb,day);
%
% fname   = name of file containing otolith profile data
% UNKname = name of series of spot scans making up otolith profile
% SRM     = name of variable in file SRM.mat to use for Standard Reference Material
%           e.g., SRM = 'nist612'
% dwell   = dwell time of the quadrupole (msec); use a scalar for constant
%           and a row vector for variable dwell times
%   adj   = optionally adjust concentrations that fall below LOD to:
%           0 (adj=0), LOD (adj=1), or make no adjustment (adj=2, default)
%   spike = optional spike removal of background-corrected signal data via Grubbs
%           Test (= 'g'), Rosner Test (='r'), or STDEV (= 1 to 4)    (default = 0)
%   drift = optionally correct for instrument drift using linear interpolation
%           (= 1) or nearest neighbor method (= 2)                   (default = 0)
%   tol   = min value of R2 required to use linear vs. nearest neighbor
%           interpolation                                         (default = 0.55)
%   verb  = verbose output of progress                               (default = 1)
%   day   = determine date of acquisition from file name             (default = 0)
%
% result = structure of results with the following fields:
%  .spot    = numeric tag identifying spot scan
%  .txt     = cell array of text for each analyte element
%  .iso     = corresponding isotope number
%  .txt_iso = combined element + isotope labels
%  .ppm     = concentration of unknown, averaged across replicates                 (ppm)
%  .LOD     = limits of detection, averaged across replicates                      (ppm)
%  .ratio   = molar ratios to internal standard, averaged across replicates (mMole/Mole)
%  .SRM     = name of SRM used for external calibration
%  .adj     = type of adjustment applied to values below LOD             (zero,LOD,none)
%  .spike   = type of spike removal applied to the time series
%  .drift   = method of drift correction
%  .tol     = tolerance for linear interpolation
%  .SRM     = name of SRM used for external calibration
%  .gDate   = cell array of Gregorian date of acquisition
%
%
%  SEE ALSO: f_plot_PS, f_cps2ppm_PT, f_cps2ppm

% -----Author:-----
% by David L. Jones, Jun-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% June-2011: added support for match, day; oto is now updated as each otolith is
%            processed; output all structure so background and LOD can be examined
%            for day effects

% -----Set defaults & check input:-----
if (nargin <  5), adj   = 2;    end % default no adjustment if PPM < LOD
if (nargin <  6), spike = 0;    end % default no spike removal
if (nargin <  7), drift = 0;    end % default no drift correction
if (nargin <  8), tol   = 0.55; end % default tolerance for linear interpolation of m
if (nargin <  9), verb  = 1;    end % default verbose output of progress
if (nargin < 10), day   = 0;    end % default don't get date from file name

% Check day:
if (day~= 0) && (day~=1)
   error('DAY must be 0 or 1!')
end

% Specify Regular Expression matching pattern:
ptn = ['^' UNKname '_\d+$']; % series of edge-to-core otolith spot scans

% Specify name of Standards:
STDname = ['cSTD_' UNKname];

% Check that file exists:
if (exist(fname,'file')~=2)
   error(['Cannot find file ' fname '!'])
end

% Make sure external standard is present:
if isempty(who('-file',fname,STDname))
   error(['File ' fname 'lacks a cell array of external standards (= ' STDname ')!']);
end

% Make sure internal standard is present:
if isempty(who('-file',fname,'IS'))
   error(['File ' fname 'needs a cell array for the internal standard (= IS)!']);
end

% Check for presence of otolith spot scans:
if isempty(who('-file',fname,'-regexp',ptn));
   error(['File ' fname ' contains no spot scans named ' UNKname '!']);
end
% -------------------------------------

% Load variables:
load(fname,STDname,'IS');  % standards
load(fname,'-regexp',ptn); % otolith spot scans

% Get list of otolith spot scans:
S.txt = who('-file',fname,'-regexp',ptn);

% -----Extract numeric spot tag:-----
S.cell = regexp(S.txt,'(?<=PS\d+_)\d+', 'match'); % spot number
nCell  = size(S.cell,1); % get # matchingcells
S.num  = nan(nCell,1);   % preallocate

for z = 1:nCell
   S.num(z) = str2num(cell2mat(S.cell{z})); % convert cell string to numbers
end

uS = unique(S.num); % list of unique spots in this file
nS = numel(uS);     % # unique spots in this file
% --------------------------------------

% Show file being processed:
if (verb>0)
   fprintf('\nProcessing file %s...\n',fname);
end

% Process each spot separately:
for j = 1:nS
   % Show spot being processed:
   if (verb>0)
      fprintf('     %s\n',S.txt{j});
   end
   
   % Process this spot:
   spot = f_cps2ppm(eval(S.txt{j}),eval(STDname),SRM,IS,dwell,adj,spike,drift,tol,0);
   
   % Get date of acquisition:
   switch day
      case 0 % Use date from f_cps2ppm (from original *.csv if available):
         acqDate = spot.gDate{:};
      case 1 % Determine date from filename:
         dateStr = cell2mat(regexp(fname,'\d{6}(?=.mat)', 'match'));
         yyVar   = str2double(dateStr(1:2))+2000;
         mmVar   = str2double(dateStr(3:4));
         ddVar   = str2double(dateStr(5:6));
         acqDate = f_gregorian(f_julian([yyVar mmVar ddVar]),yyVar);
   end
   
   % Collect data for this spot:
   if (j==1)
      % Initialize with first spot:
      all.spot    = S.num(j);
      all.txt     = spot.txt;
      all.iso     = spot.iso;
      all.txt_iso = spot.txt_iso;
      all.ppm     = spot.ppm;
      all.LOD     = spot.LOD;
      all.ratio   = spot.ratio;
      all.gDate   = acqDate;
   else
      all.spot  = [all.spot;  S.num(j)];
      all.ppm   = [all.ppm;   spot.ppm];
      all.LOD   = [all.LOD;   spot.LOD];
      all.ratio = [all.ratio; spot.ratio];
      all.gDate = [all.gDate; acqDate];
      
      % Make sure samples from all spots include the same isotopes:
      if (~isequal(all.txt_iso,spot.txt_iso))
         error([S.txt{j} ' does not include the same isotopes as the previous sample!']);
      end
   end
end
fprintf('Done!\n');

% Rename for output:
result = all;
