function f_graphviz_neato(fname,PDF)
% - create an undirected graph using Graphviz
% 
% USAGE: f_graphviz_neato('fname',PDF);
% 
% fname = name of '*.dot' file to process
% PDF   = generate PDF file in Preview.app              (default = 1)
% 
% SEE ALSO: f_graphviz_mds, f_graphviz_mst

% -----Notes:----
% This function requires the Graphviz program, which is available from:
% http://www.graphviz.org
% 
% Depending on your platform, you may also need to: (1) edit the
% 'f_graphviz_neato' function to specify the correct location of 'neato' defined
% in the 'pname' variable; and (2) comment out the section that calls the
% external 'Preview.app' if you are not on OS X. Windows users could call Adobe
% Distiller instead and Linux users could call Ghostscript.

% -----Author:-----
% by David L. Jones, Jun-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), PDF = 1; end % generate PDF by default

% Check operating system:
if (ismac<1)
   error('This function currently requires MacOS X!')
end

% Check that file 'fname' is in your path:
if (exist(fname,'file')~=2)
   error(['Can''t find file ' fname '!'])
end

% Set program path & file name:
pname = '/usr/local/bin/neato';
% -------------------------------------

% Set name of output file:
cmd      = '(?<name>.+?)\.dot';
S        = regexp(fname,cmd,'names');
fnameDOT = [S.name '.dot']; % DOT file to import
fnameEPS = [S.name '.eps']; % Postscript file to export
clear cmd S:

% Call NEATO to process DOT file:
eval(['!' pname ' -Teps ' fnameDOT ' -o ' fnameEPS]);

% Optionally display PDF in Preview.app:
if (PDF>0)
   eval(['!open -a Preview ' fnameEPS]);
end
