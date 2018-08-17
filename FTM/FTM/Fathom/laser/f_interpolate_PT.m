function result = f_interpolate_PT(tra,vNames)
% - interpolate otolith profile transects so they are equal lengths
%
% USAGE: result = f_interpolate_PT(tra,vNames);
%
% tra    = cell array of otolith profile transects created by f_cps2ppm_PT or
%          f_extract_PT (e.g., tra = {sub_1 sub_2 sub_3})
% vNames = cell array of output variable names
%          (e.g., tra = {'sub_1i' 'sub_2i' 'sub_3i'})
% 
% results = optional structure of results with the following fields:
%  .meanPT = mean otolith profile transect
%  .SD     = corresponding standard deviation
% 
% SEE ALSO: f_compare_PT, f_extract_PT, f_cps2ppm_PT

% -----Notes:-----
% This function makes several otolith profile transects compatible for making
% direct comparisons (e.g., calculating deviations, generating mean profiles, or
% determining r^2) by interpolating all profiles so they all have the same
% number of data points as the longest profile. Note interpolated versions of
% the input variables are exported directly to the MATLAB workspace according to
% the names provided in by the VNAMES variable. Additional output is also
% optionally returned when the output argument RESULTS is specified. 

% -----Author:-----
% by David L. Jones, Aug-2013
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% Sep-2014: updated documentation

% -----Set defaults and check input:-----
if (~iscell(tra) || numel(tra)<2)
   error('TRA must be a cell array of multiple otolith profile transects!');
end
if (~iscell(vNames))
   error('vName must be a cell array specifying output variable names');
end
if (~isequal(numel(tra),numel(vNames)))
   error('VNAME and TRA are not of compatible sizes!');
end
% ---------------------------------------

% Get length of otolith transects:
nTra = numel(tra);
n    = zeros(nTra,1); % preallocate
c    = zeros(nTra,1);
for i=1:nTra
   n(i) = numel(tra{i}.pos);  % get # rows of each transect
   c(i) = size(tra{i}.ppm,2); % get # cols of each transect
end

% Get # cols:
nc = f_unique(c);
if (numel(nc)>1)
   error('All members of TRA must have the same # columns in the *.ppm matrix!');
end

% Get # rows of longest transect:
nr = max(n);

% Interpolate values so all transects are equal length:
for i=1:nTra
   nc             = size(tra{i}.ppm,2); % get # columns
   tra{i}.pos_i   = linspace(tra{i}.pos(1),tra{i}.pos(end),nr);
   tra{i}.ppm_i   = nan(nr,nc);         % preallocate
   tra{i}.LOD_i   = nan(nr,nc);
   tra{i}.ratio_i = nan(nr,nc);
   for j = 1:nc
      tra{i}.ppm_i(:,j)   = interp1(tra{i}.pos,tra{i}.ppm(:,j),tra{i}.pos_i);
      tra{i}.LOD_i(:,j)   = interp1(tra{i}.pos,tra{i}.LOD(:,j),tra{i}.pos_i);
      tra{i}.ratio_i(:,j) = interp1(tra{i}.pos,tra{i}.ratio(:,j),tra{i}.pos_i);
   end
end

% -----Calculate mean otolith profile transect:-----
ppm_mean   = nan(nr,nc); % preallocate
LOD_mean   = nan(nr,nc);
ratio_mean = nan(nr,nc);
ppm_SD     = nan(nr,nc);
LOD_SD     = nan(nr,nc);
ratio_SD   = nan(nr,nc);

for i=1:nc % repeat for each mass
   ppm_temp   = nan(nr,nTra); % preallocate
   LOD_temp   = nan(nr,nTra);
   ratio_temp = nan(nr,nTra);
   for j=1:nTra % repeat for each transect
      ppm_temp(:,j)  = tra{j}.ppm_i(:,i);
      LOD_temp(:,j)  = tra{j}.LOD_i(:,i);
      ratio_temp(:,j) = tra{j}.ratio_i(:,i);
   end
   
   % Calculate means (and STDVs) of this mass across all transects:
   ppm_mean(:,i)   = mean(ppm_temp,2);
   LOD_mean(:,i)   = mean(LOD_temp,2);
   ratio_mean(:,i) = mean(ratio_temp,2);
   
   ppm_SD(:,i)     = std(ppm_temp,0,2);
   LOD_SD(:,i)     = std(LOD_temp,0,2);
   ratio_SD(:,i)   = std(ratio_temp,0,2);
end
% --------------------------------------------------

% Export interpolated data directly MATLAB workspace:
for i=1:nTra
   % 'assignin' can't work with structures, so:
   assignin('base','assignin_temp_txt',tra{i}.txt);
   assignin('base','assignin_temp_iso',tra{i}.iso);
   assignin('base','assignin_temp_txt_iso',tra{i}.txt_iso);
   assignin('base','assignin_temp_pos',tra{i}.pos_i);
   assignin('base','assignin_temp_ppm',tra{i}.ppm_i);
   assignin('base','assignin_temp_LOD',tra{i}.LOD_i);
   assignin('base','assignin_temp_ratio',tra{i}.ratio_i);
   
   evalin('base',[vNames{i} '.txt     = assignin_temp_txt;']);
   evalin('base',[vNames{i} '.iso     = assignin_temp_iso;']);
   evalin('base',[vNames{i} '.txt_iso = assignin_temp_txt_iso;']);
   evalin('base',[vNames{i} '.pos     = assignin_temp_pos;']);
   evalin('base',[vNames{i} '.ppm     = assignin_temp_ppm;']);
   evalin('base',[vNames{i} '.LOD     = assignin_temp_LOD;']);
   evalin('base',[vNames{i} '.ratio   = assignin_temp_ratio;']);
   evalin('base','clear assignin_temp*')
end

% Wrap results up into a structure:
result.ppm_mean   = ppm_mean;
result.ppm_SD     = ppm_SD;
result.LOD_mean   = LOD_mean;
result.LOD_SD     = LOD_SD;
result.ratio_mean = ratio_mean;
result.ratio_SD   = ratio_SD;
