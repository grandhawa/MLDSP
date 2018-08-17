% f_import - imports comma-delimited file having row and col labels
%
% Usage: f_import
% 
% ----------Notes:------------------------------------------
% This script allows you to select a comma-delimited ASCII
% file via a GUI for import into a Matlab matrix. The row &
% column labels are converted to cell arrays. The REQUIRED
% format of the ASCII file is shown in the example.
%
% -----Example ASCII file for import:-----
% stn:1,stn:2,stn:3,stn:4,stn:5
% sp_A,11,22,33,44,55
% sp_B,0.2,0.6,0.7,0.1,0.0
% sp_C,122,66,902,7,66
% sp_D,0.9,33,0.04,5,99
% ------------------------------------------
%
% Here we have abundances of 4 species (rows) from 5 stations (cols)
%
% SEE ALSO: f_export

% -----Author:-----
% by David L. Jones, Feb-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% ----- Dependencies: ---------------------------------
% requires tblread.m from the TMW's Statistics Toolbox
% -----------------------------------------------------

% 19-July-2001: added label conversion to cell arrays

[fileVar,pathVar] = uigetfile('*.d??', 'Specify file to import');
nameVar = fileVar(1,1:(size(fileVar,2) - 4));
colVar = [nameVar '_col'];
rowVar = [nameVar '_row'];
% use 2 quotes to specify a quote in a string:
commandVar = ['[' nameVar ',' colVar ',' rowVar '] = tblread([''' pathVar fileVar '''],''comma'');'];

eval(commandVar); % this is what is returned

% convert row & col labels to cell arrays:
commandVar = [colVar ' = cellstr(' colVar ');'];
eval(commandVar);

commandVar = [rowVar ' = cellstr(' rowVar ');'];
eval(commandVar);

clear fileVar pathVar rowVar colVar commandVar nameVar; % cleanup...not needed if this was a function
whos; % list variables in workspace
