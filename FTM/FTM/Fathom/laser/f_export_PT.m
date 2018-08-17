function f_export_PT(fname)
% - export otolith profile data to a comma-separated-values (CSV) text file
%
% USAGE: f_export_PT('fname')
%
% fname = cell array of file(s) to process
%
% SEE ALSO: f_cps2ppm_PT, f_plot_PT, f_export_cps

% -----Author:-----
% by David L. Jones, Oct-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Dec-2011: now exports all data types (PPM, RATIO, and LOD)
% Apr-2011: no longer exports RATIO; now exports raw CPS when available

% -----Set defaults & check input:-----
if (nargin <  2), verb  = 1; end % default verbose output of progress

ptn = '^tra_\w+$'; % match any variable beginning with "tra_"

% Check that FNAME is a cell array:
if ~iscell(fname), error('FNAME must be a cell array'); end

% Check that files in FNAME exist and contain the proper variables:
nFiles = numel(fname); % get # files
for i = 1:nFiles
   % Check that file exists:
   if (exist(fname{i},'file')~=2)
      error(['Cannot find file ' fname(i) '!'])
   end
   
   % Check for presence of otolith profile transects:
   if isempty(who('-file',fname{i},'-regexp',ptn));
      error(['File ' fname{i} ' contains no profile transects!']);
   end
end
% -------------------------------------

% Process each file separately:
for i = 1:nFiles
   % Load variables:
   load(fname{i},'-regexp',ptn); % otolith profile transects
   
   % Get list of otolith profile transects:
   PT.txt = who('-file',fname{i},'-regexp',ptn);
   nPT    = size(PT.txt,1); % get # profile transects in this file
   
   % Show file being processed:
   if (verb>0)
      fprintf('\nProcessing file %s...\n',fname{i});
   end
   
   % Process each profile transect separately:
   for j = 1:nPT
      X = eval(PT.txt{j});
      
      % Setup export filename extensions:
      out_PPM = [PT.txt{j} '_ppm.csv'];
      out_LOD = [PT.txt{j} '_lod.csv'];
      out_CPS = [PT.txt{j} '_cps.csv'];
      
      % Don't overwrite existing file:
      if exist(out_PPM,'file'), error(['File ' out_PPM ' already exists!']); end
      if exist(out_LOD,'file'), error(['File ' out_LOD ' already exists!']); end
      if exist(out_CPS,'file'), error(['File ' out_CPS ' already exists!']); end
      
      % -----Export data:-----
      % PPM:
      f_exportR([X.pos X.ppm],cellstr(num2str([1:numel(X.s)]')),['pos' X.txt_iso],out_PPM);
      
      % LOD:
      f_exportR([X.pos X.LOD],cellstr(num2str([1:numel(X.s)]')),['pos' X.txt_iso],out_LOD);
      
      % RAW CPS:
      if ~isnan(X.cps)
         f_exportR([X.pos X.cps ],cellstr(num2str([1:numel(X.s)]')),['pos' X.txt_iso],out_CPS);
      end
   end
   
   % Clear variables associated with this file:
   clear PT;
   clear('-regexp',ptn);
end
