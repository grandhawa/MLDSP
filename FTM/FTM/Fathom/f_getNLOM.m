function f_getNLOM(x,y,region)
% - download Navy Layered Ocean Model SSH Nowcast imagery
%
% USAGE: f_getNLOM(x,y,region)
%
% x      = first day of imagery to download (e.g., [yyyy mm dd])
% y      = last day
% region = 1:IAS, 2:GOM

% -----NOTES:-----
% This function downloads the images to the current working directory, so be
% sure to change to the appproriate directory before running.

% -----Author:-----
% by David L. Jones, Jun-2006
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Sep-2006: added support for GOM

% Set up default region (Inter-American Seas, Nowcast, SSH):
urlVar    = ('http://www7320.nrlssc.navy.mil/global_nlom32/navo/');
if (region==1)
	headerVar = ('IAS/TOPEX+ERS2+GFO+MOD_SSH_IAS_');
else
	headerVar = ('GOM/TOPEX+ERS2+GFO+MOD_SSH_GOM_');
end


% Create vector of date numbers:
nDates = (datenum(x):1:datenum(y))';

% Create vector of date strings
sDates = datestr(nDates,29);

% Remove the hyphens ('-') from date strings
sDates(:,5) = '';
sDates(:,7) = '';

% Get number of dates:
noDates = size(sDates,1);

% Download images to present working directory:
for i = 1:noDates
   fprintf(['Downloading file ' num2str(i) ' of ' num2str(noDates) '\n'])
   urlwrite([urlVar headerVar sDates(i,:) '.001.gif'],[sDates(i,:) '.gif']);
end

fprintf(['\n\n Downloading complete...' '\n\n'])
