function [result,all] = f_cps2ppm_SPOT(fname,SRM,IS,dwell,adj,spike,drift,tol,verb,match,day)
% - batch process parsed LA-ICP-MS spot data
%
% USAGE: [result,all] = f_cps2ppm_SPOT(fname,'SRM',IS,dwell,adj,spike,drift,tol,verb,match,day);
%
% fname   = cell array of file(s) to process
% SRM     = name of variable in file SRM.mat to use for Standard Reference Material
%           e.g., SRM = 'nist612'
% IS      = structure of Internal Standard with the following fields:
%           .txt  = cell array of text indicating element
%           .ppm  = corresponding concentration in Unknown
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
%   match = regular expression pattern to match names of otolith spot scans
%           1: matches o1140_A (replicates A-C)                      (default = 1)
%           2: matches o1140_A (all replicates)
%           3: matches oto10A1 
%           4: matches f1140_A (replicates A-C) 
%           5: matches s1140_A (replicates A-C) 
%   day   = determine date of acquisition from file name 
%
% result = structure of results with the following fields:
%  .oto     = numeric tag identifying individual otoliths
%  .txt     = cell array of text for each analyte element
%  .iso     = corresponding isotope number
%  .txt_iso = combined element + isotope labels
%  .nRep    = # within-otolith replicate spot samples used to calclate averages
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
% all = structure of results (NOT averaged acrossed otoliths):
%  .oto     = numeric tag identifying individual otoliths
%  .txt     = cell array of text for each analyte element
%  .iso     = corresponding isotope number
%  .txt_iso = combined element + isotope labels
%  .ppm     = concentration of unknown                                     (ppm)
%  .LOD     = limits of detection                                          (ppm)
%  .ratio   = molar ratios to internal standard                     (mMole/Mole)
%  .bg      = mean background counts of UNK                                (cps)
%  .gDate   = cell array of Gregorian date of acquisition
%
%  SEE ALSO: f_cps2ppm, f_cps2ppm_MACS, f_importXL, f_cpsParse

% -----Notes:-----
% This function is used to process all the LA-ICP-MS spot sample data
% within the file(s) specified by FNAME. Each file should have
% already been processed with f_importXL and f_cpsParse and must contain the
% following variables:
%
% cSTD, fSTD, and/or sSTD = cell array of multiple external standards, each
%                           created by f_cpsParse.
% 
% One of several types of structures created by f_cpsParse adhering to the
% following naming convention:
%
% o####_A:  1) begins with the letter 'o' (= otolith)
%           2) followed by a variable length numeric tag unique to each fish
%           3) followed by an underscore and letter (A-C) indicating the
%              replicate  spot sample taken from the otolith of this fish
%           Example: o1120_B refers to the variable containing the second
%           replicate of otolith spot sample data taken from fish # 1120
%
% oto####A1: 1) begins with the letters 'oto' (= otolith)
%            2) followed by a variable length numeric tag unique to each fish
%            3) followed a letter indicating right vs. left otolith (A or B)
%            4) followed by a number indicating the replicate spot sample
%               taken from the otolith of this fish
%            Example: oto10B2 refers to the variable containing the second
%            replicate of otolith spot sample data taken from fish # 10
% 
% f####_A:   similar to 'o####_A', but used to match fin rays rather than
%            otoliths; requires the corresponding external standards 
%            to  be named 'fCTD' and an appropriate internal standard (IS)
% 
% s####_A:   similar to 'o####_A', but used to match fin spines rather
%            than otoliths; requires the corresponding external standards
%            to  be named 'sCTD' and an appropriate internal standard (IS)
%
% By default, the date of acquisition (gDate) is determined from the UNK
% structure passed to the f_cps2ppm function, which may have originated
% from the original *.csv file imported by the f_importXL function.
% However, if the acquistion date is missing (e.g., gDate = {'NaN'}), it can be
% determined from the filename, if they are named according to the 'yymmdd.mat'
% format, by setting DAY=1.
%
% Setting VERB=1 assists in tracking down which files and/or variables
% produce errors in f_cps2ppm when, say, no background is provided or the
% background and signal portions overlap, etc.
%
% Replicate spot scans from each sample are first checked for the presence of
% outliers (based on PPM), then within-sample averages of PPM, LOD, and
% RATIO are calculated for those replicates not identified as outliers.

% -----Author:-----
% by David L. Jones, May-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jun-2011: added support for match, day; oto is now updated as each otolith is
%           processed; output all structure so background and LOD can be examined
%           for day effects
% Apr-2012: updated error messages
% May-2012: added support for samples with different internal standards, so
%           now you must specify 'IS' on command line not in data file; now
%           supports files cell array so external standards named cSTD,
%           fSTD, and/or sSTD.
% Jun-2012: don't clear xSTD or IS after each file is processed
% Aug-2012: now clears xSTD after each file is processed; added 'all.fname'
%           field in output; fixed "file "exist" error checking

% -----Set defaults & check input:-----
if (nargin <  5), adj   = 2;    end % default no adjustment if PPM < LOD
if (nargin <  6), spike = 0;    end % default no spike removal
if (nargin <  7), drift = 0;    end % default no drift correction
if (nargin <  8), tol   = 0.55; end % default tolerance for linear interpolation of m
if (nargin <  9), verb  = 1;    end % default verbose output of progress
if (nargin < 10), match = 1;    end % default regular expression pattern
if (nargin < 11), day   = 0;    end % default don't get date from file name

% Specify Regular Expression matching pattern:
switch match
   case 1
      ptn  = '^o\d+_[A-C]$';  % Dave's otolith naming convention (only reps A-C)
      xSTD = 'cSTD';          % name of cell array of external standards
   case 2
      ptn  = '^o\d+_[A-Z]+$'; % Dave's otolith naming convention (all replicates)
      xSTD = 'cSTD';          % name of cell array of external standards
   case 3
      ptn  = '^oto\d+[A-B]\d$'; % Holly's otolith naming convention
      xSTD = 'cSTD';            % name of cell array of external standards
   case 4
      ptn = '^f\d+_[A-C]$';   % Ori's fin ray naming convention (only reps A-C)
      xSTD = 'fSTD';          % name of cell array of external standards for FINS
   case 5
      ptn = '^s\d+_[A-C]$';   % Ori's fin spine naming convention (only reps A-C)
      xSTD = 'sSTD';          % name of cell array of external standards for SPINES
   otherwise
      error('Unknown regex pattern specified by PTN!')
end

% Check day:
if (day~= 0) && (day~=1)
   error('DAY must be 0 or 1!')
end
% -------------------------------------

% Check that FNAME is a cell array:
if ~iscell(fname), error('FNAME must be a cell array'); end

% Check that files in FNAME exist and contain the proper variables:
nFiles = numel(fname); % get # files
for i = 1:nFiles
   % Check that file exists:
   if (exist(fname{i},'file')~=2)
      error(['Cannot find file ' fname{i} '!'])
   end
   
   % Make sure external standard is present:
   if isempty(who('-file',fname{i},xSTD))
      error(['File ' fname{i} ' needs a cell array for the external standards (= ' xSTD ')!']);
   end
        
   % Check for presence of otolith spot scans:
   if isempty(who('-file',fname{i},'-regexp',ptn));
      error(['File ' fname{i} ' contains no otolith spot scans!']);
   end
end
% -------------------------------------

% Initialize variables:
oto     = [];
txt     = [];
iso     = [];
txt_iso = [];
nRep    = [];
ppm     = [];
LOD     = [];
ratio   = [];
gDate   = [];
cnt     = 0;

% Process each file separately:
for i = 1:nFiles
   
   % Show file being processed:
   if (verb>0)
      fprintf('\nProcessing file %s...\n',fname{i});
   end
   
   % Load variables:
   load(fname{i},xSTD);          % standards
   load(fname{i},'-regexp',ptn); % spot scans
   
   % Get list of spot scans:
   S.txt = who('-file',fname{i},'-regexp',ptn);
   
   % -----Extract numeric sample tag:-----
   switch match
      case 1
         S.cell = regexp(S.txt,'(?<=o)\d+(?=_[A-C])', 'match');  % reps A-C
      case 2
         S.cell = regexp(S.txt,'(?<=o)\d+(?=_[A-Z]+)', 'match'); % all replicates
      case 3
         S.cell = regexp(S.txt,'(?<=oto)\d+(?=[A-B]\d)', 'match');
      case 4
         S.cell = regexp(S.txt,'(?<=f)\d+(?=_[A-C])', 'match');  % reps A-C
      case 5
         S.cell = regexp(S.txt,'(?<=s)\d+(?=_[A-C])', 'match');  % reps A-C
   end
   nCell = size(S.cell,1); % get # matching cells
   S.num = nan(nCell,1);   % preallocate
   for z = 1:nCell
      S.num(z) = str2num(cell2mat(S.cell{z})); % convert cell string to numbers
   end
   uS = unique(S.num); % list of unique otoliths in this file
   nS = numel(uS);     % # unique otoliths in this file
   % --------------------------------------
   
   % Process each otolith separately:
   for j = 1:nS
      idxS      = find(S.num == uS(j)); % get index to spot scans matching this otolith
      rep       = numel(idxS);          % get # replicates for this otolith
      oto_ppm   = [];                   % initialize
      oto_LOD   = [];
      oto_ratio = [];
      
      % Collect list of otoliths across files:
      oto = [oto; uS(j)];
      
      % Process each spot scan separately:
      for k = 1:rep
         
         % Show spot being processed:
         if (verb>0)
            fprintf('     %s\n',S.txt{idxS(k)});
         end
                          
         % Process this replicate:
         spot = f_cps2ppm(eval(S.txt{idxS(k)}),eval(xSTD),SRM,IS,dwell,adj,spike,drift,tol,0);
                          
         % Get date of acquisition:
         switch day
            case 0 % Use date from f_cps2ppm (from original *.csv if available):
               acqDate = spot.gDate{:};
            case 1 % Determine date from filename:
               dateStr = cell2mat(regexp(fname{i},'\d{6}(?=.mat)', 'match'));
               yyVar   = str2double(dateStr(1:2))+2000;
               mmVar   = str2double(dateStr(3:4));
               ddVar   = str2double(dateStr(5:6));
               acqDate = f_gregorian(f_julian([yyVar mmVar ddVar]),yyVar);
         end
         
         % Add variable labels only once:
         if isempty(txt)
            txt     = spot.txt;
            iso     = spot.iso;
            txt_iso = spot.txt_iso;
         else
            % Make sure otoliths from all files include the same isotopes:
            if (~isequal(txt_iso,spot.txt_iso))
               error([S.txt{idxS(k)} ' does not include the same isotopes as the previous sample!']);
            end
         end
         
         % Collect spot data for this otolith:
         oto_ppm   = [oto_ppm;   spot.ppm];
         oto_LOD   = [oto_LOD;   spot.LOD];
         oto_ratio = [oto_ratio; spot.ratio];
         
         % -----Optionally collect all data (no averaging):-----
         if (nargout>1)
            % Record the file it came from:
            cnt            = cnt + 1;
            all.fname{cnt} = fname{i};
         
            if (i==1 && j==1 && k==1)
               % Initialize with first spot:
               all.oto     = uS(j);
               all.txt     = spot.txt;
               all.iso     = spot.iso;
               all.txt_iso = spot.txt_iso;
               all.ppm     = spot.ppm;
               all.LOD     = spot.LOD;
               all.ratio   = spot.ratio;
               all.bg      = spot.bg;
               all.gDate   = acqDate;
            else
               all.oto     = [all.oto;   uS(j)];
               all.ppm     = [all.ppm;   spot.ppm];
               all.LOD     = [all.LOD;   spot.LOD];
               all.ratio   = [all.ratio; spot.ratio];
               all.bg      = [all.bg;    spot.bg];
               all.gDate   = [all.gDate; acqDate];
            end
         end
         % -----------------------------------------------------
      end
      
      % Check for outliers:
      if (rep>2)
         out = f_grpOutlier(oto_ppm,ones(rep,1),10); % 0 = outlier
      else
         out = ones(rep,1);                          % no outliers
      end
      
      % Average values across replicates (excluding outliers):
      ppm   = [ppm;   mean(oto_ppm(out,:),1)];
      LOD   = [LOD;   mean(oto_LOD(out,:),1)];
      ratio = [ratio; mean(oto_ratio(out,:),1)];
      nRep  = [nRep; sum(out>0)]; % get # replicates used for averages
      
      % Add date of acquisition:
      gDate = [gDate; acqDate];
   end
   
   % Clear variables associated with this file:
   clear S;
   clear('-regexp',ptn);
   eval(['clear ' xSTD])
end

% Wrap results up into a structure:
result.oto     = oto;        % numeric tag identifying individual otoliths
result.txt     = txt;        % cell array of text for each analyte element
result.iso     = iso;        % corresponding isotope number
result.txt_iso = txt_iso;    % combined element + isotope labels
result.nRep    = nRep;       % within-otolith replicate spot samples used to calclate averages
result.ppm     = ppm;        % concentration of unknown, averaged across replicates (ppm)
result.LOD     = LOD;        % limits of detection, averaged across replicates (ppm)
result.ratio   = ratio;      % molar ratios to internal standard, averaged across replicates (mMole/Mole)
result.SRM     = spot.SRM;   % name of SRM used for external calibration
result.adj     = spot.adj;   % type of adjustment applied to values below LOD
result.spike   = spot.spike; % type of spike removal applied to the time series
result.drift   = spot.drift; % type of drift applied
result.tol     = spot.tol;   % tolerance for linear interpolation
result.gDate   = gDate;      % date of acquisition

