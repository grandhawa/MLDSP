function f_export_cps(X,fname)
% - export LA-ICP-MS transient signal cps data
%
% USAGE: f_export_cps(X,'fname')
%
% X     = structure of measured SAMPLE signal previously created by f_cpsParse
% fname = name of destination file (omit *.csv extension)
%
% SEE ALSO: f_cpsParse, f_export_PT

% -----Author:-----
% by David L. Jones, Oct-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
% Don't overwrite existing file:
if exist([fname '.csv'],'file')
   error('Destination file already exists!')
end
% -------------------------------------

% Add combined element + isotope labels:
txt_iso = cellstr(strcat(X.txt',regexprep( cellstr(num2str(X.iso')),'\s', '') ))';

% Export SIGNAL data as CSV:
f_exportR([X.s X.cps],cellstr(num2str([1:numel(X.s)]')),...
   ['sec' txt_iso],[fname '.csv']);

% Export BACKGROUND data as CSV:
f_exportR([X.bg.s X.bg.cps],cellstr(num2str([1:numel(X.bg.s)]')),...
   ['sec' txt_iso],[fname '_background.csv']);
