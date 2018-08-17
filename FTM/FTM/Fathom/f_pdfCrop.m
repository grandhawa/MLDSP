function f_pdfCrop(input,output)
% - crop a PDF file to its bounding box
%
% USAGE: f_pdfCrop('input.pdf','output.pdf')
% 
% SEE ALSO: f_pdf, f_pdfDistill, f_pdfLatex, f_pdfMerge, f_pdfSplit

% -----Notes:-----
% 'pdfcrop.pl' is an executable PERL script located in '~/bin/' which calls
% Ghostscript to determine the bounding box size of each page in a PDF, then
% crops the PDF accordingly. However, the PERL script cannot be called directly
% using the Matlab 'system' or bang ('!') command because when PERL calls
% Ghostscript the script fails. Since the script can be called successfully
% outside the Matlab environment (e.g., from Terminal) it is likely the Matlab
% runtime environment does not provide the sufficient paths, etc. for PERL to
% call Ghostscript. So, this function is used to call the PERL script from a
% Terminal window via an AppleScript.

% http://apple.stackexchange.com/questions/24431/osascript-error
% 
% On 64-bit sytems you may need to specify the architecture to run the
% AppleScript loader is 32 bit:
% Error loading....Library/ScriptingAdditions/YouHelper.osax/Contents/MacOS/YouHelper

% -----Author:-----
% by David L. Jones, Dec-2008
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Feb-2014: added fix for 'arch -i386'

% -----Check input & set defaults:-----
if (nargin < 2), output = input; end % default use same name for cropped version

% Add full directory name:
input  = [pwd '/' input];
output = [pwd '/' output];

% Build an AppleScript:
% cmd_1 = ['osascript -e ''tell application "Terminal"'' '];
cmd_1 = ['arch -i386 osascript -e ''tell application "Terminal"'' ']; % force 32 bit
cmd_2 = ['-e ''activate'' '];
cmd_3 = ['-e ''do script "pdfcrop ' input ' ' output '"'' '];
cmd_4 = ['-e ''end tell'' '];

% Invoke Applescript from Matlab:
system([cmd_1 cmd_2 cmd_3 cmd_4]);



