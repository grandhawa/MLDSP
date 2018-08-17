function result = f_cps2ppm_MACS(fname,SRM,dwell,adj,spike,drift,tol,verb,day)
% - batch process parsed LA-ICP-MS spot data for MACS-3 samples
%
% USAGE: result = f_cps2ppm_MACS(fname,'SRM',dwell,adj,spike,drift,tol,verb,day);
%
% fname   = cell array of file(s) to process
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
%
% result = structure of results with the following fields:
%  .macs    = tag identifying individual spot samples
%  .txt     = cell array of text for each analyte element
%  .iso     = corresponding isotope number
%  .txt_iso = combined element + isotope labels
%  .ppm     = concentration of unknown                                     (ppm)
%  .LOD     = limits of detection                                          (ppm)
%  .ratio   = molar ratios to internal standard                     (mMole/Mole)
%  .bg      = mean background counts of UNK                                (cps)
%  .gDate   = cell array of Gregorian date of acquisition
%
% SEE ALSO: f_cps2ppm, f_cps2ppm_SPOT

% -----Notes:-----
% This function is used to process all the LA-ICP-MS MACS-3 spot sample data
% within the file (or files) specified by FNAME. Each file should have
% already been processed with f_importXL and f_cpsParse and requires the
% following variables:
%
% cSTD    = cell array of multiple STD's, each created by f_cpsParse
% macs_A  = structure created by f_cpsParse adhering to the following
%           naming convention:
%           1) begins with the letters 'macs' (= MACS-3 carbonate standard)
%           2) followed by an underscore and letter indicating the unique
%              spot sample taken during this day's run
%           Example: macs_C refers to the variable containing the third
%           MACS-3 spot sample data taken during this day's run
%
% Setting verb=1 assists in tracking down which files and/or variables
% produce errors in f_cps2ppm when, say, no background is provided or the
% background and signal portions overlap, etc.

% -----Author:-----
% by David L. Jones, June-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jan-2015: corrected in 'exist' error message; removed '.S' field for
%           normalized sensitivity since f_cps2ppm no longer returns this either

% -----Set defaults & check input:-----
if (nargin <  4), adj   =  2;    end % default no adjustment if PPM < LOD
if (nargin <  5), spike =  0;    end % default no spike removal
if (nargin <  6), drift =  0;    end % default no drift correction
if (nargin <  7), tol   =  0.55; end % default tolerance for linear interpolation of m
if (nargin <  8), verb  =  1;    end % default verbose output of progress
if (nargin <  9), day   = 0;     end % default don't get date from file name

% Specify Regular Expression matching pattern:
ptn = '^macs_[A-Z]+$';

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
   if isempty(who('-file',fname{i},'cSTD'))
      error(['File ' fname{i} ' needs a cell array for the external standards (= cSTD)!']);
   end
   
   % Check for presence of MACS-3 spot scans:
   if isempty(who('-file',fname{i},'-regexp',ptn));
      error(['File ' fname{i} 'contains no MACS-3 spot scans!']);
   end
end
% -------------------------------------

% Setup Internal Standard (IS):
IS_macs.txt = {'Ca'};
IS_macs.ppm = 376900; % MACS-3

% Process each file separately:
for i = 1:nFiles
   % Load variables:
   load(fname{i},'cSTD');        % standards
   load(fname{i},'-regexp',ptn); % MACS-3 spot scans
   
   % Get list of MACS spot scans:
   S.txt{i} = who('-file',fname{i},'-regexp',ptn);
   nS       = numel(S.txt{i}); % # spot scans in this file

   % Show file being processed:
   if (verb>0)
      fprintf('\nProcessing file %s...\n',fname{i});
   end
   
   % Process each MACS-3 spot separately:
   for j = 1:nS
      
      % Show spot being processed:
      if (verb>0)
         fprintf('     %s\n',S.txt{i}{j});
      end
      
      % Process this spot:
      spot = f_cps2ppm(eval(S.txt{i}{j}),cSTD,SRM,IS_macs,dwell,adj,spike,drift,tol,0);
      
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
      
      % Collect MACS data for this spot:
      if (i==1 && j==1)
         % Initialize with first spot:
         all.macs    = S.txt{i}{j};
         all.txt     = spot.txt;
         all.iso     = spot.iso;
         all.txt_iso = spot.txt_iso;
         all.ppm     = spot.ppm;
         all.LOD     = spot.LOD;
         all.ratio   = spot.ratio;
         all.bg      = spot.bg;
         all.gDate   = acqDate;
      else
         all.macs  = [all.macs;  S.txt{i}(j)];
         all.ppm   = [all.ppm;   spot.ppm];
         all.LOD   = [all.LOD;   spot.LOD];
         all.ratio = [all.ratio; spot.ratio];
         all.bg    = [all.bg;    spot.bg];
         all.gDate = [all.gDate; acqDate];
         
         % Make sure MACS-3 samples from all files include the same isotopes:
         if (~isequal(all.txt_iso,spot.txt_iso))
            error([S.txt{j} ' does not include the same isotopes as the previous sample!']);
         end
      end
   end
   % Clear variables associated with this file:
   clear cSTD S;
   clear('-regexp',ptn);
end

% Rename for output:
result = all;
