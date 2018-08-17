function result = f_corr_PT(X,Y,verb)
% - linear correlation between masses for replicate otolith profile transects
% 
% USAGE result = f_corr_PT(X,Y)
% 
% X,Y  = strctures created by the f_cps2ppm_PT function
% verb  = optionally send results to display                 (default = 1)
% 
% result = structure with the following fields:
%  .txt_iso = combined element + isotope labels for each mass
%  .corr    = linear correlation between corresponding masses
% 
% SEE ALSO: f_cps2ppm_PT, f_cps2ppm, f_corr

% -----Notes:-----
% This function is used to compare 2 core-to-edge otolith profile transects
% run repeatedly over the same surface. This is used to determine whether
% a consistent signal is obtained from multiple passes over the same
% otolith material.

% -----Author:-----
% by David L. Jones, Dec-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 3), verb  =  1; end % send output to display by default

% Make sure both inputs measure the same masses:
if ~isequal(X.txt_iso,Y.txt_iso)
   error('X and Y don''t contain data for the same masses');
end
% -------------------------------------

% Make sure they are the same length (assuming both had similar start locations):
nr  = min( [size(X.ppm,1) size(Y.ppm,1)] );
X_ppm = X.ppm(1:nr,:);
Y_ppm = Y.ppm(1:nr,:);

% Calculate linear correlations among corresponding elements:
nc = size(X_ppm,2); % get # of masses
corr = nan(1,nc);   % preallocate
for i = 1:nc
   corr(i) = f_corr(X_ppm(:,i),Y_ppm(:,i),0,0);
end

% Wrap results up into a structure:
result.txt_iso = X.txt_iso;
result.corr    = corr;

% -----Send output to display:-----
if (verb>0)
   table.txt_iso = result.txt_iso;
   table.corr    = result.corr;
     
   headerCell = {'Analyte','R'}; % Set up column labels
   resultCell = f_struct2flat(table);       % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat
   
   fprintf('\n==============================================================================\n');
   disp(resultCell);
   fprintf('------------------------------------------------------------------------------\n');
   disp('R is the linear correlation between corresponding analytes')
end
% ---------------------------------
