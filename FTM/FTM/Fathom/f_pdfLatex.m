function f_pdfLatex(fname,verb)
% - typeset a LaTeX file
% 
% USAGE = f_pdfLatex('fname',verb)
% 
% fname   = name of input file (e.g., = 'layout.tex')
% verb    = display merged PDF in Preview.app         (default = 0)
% 
% SEE ALSO: f_pdf, f_pdfCrop, f_pdfDistill, f_pdfMerge, f_pdfSplit


% -----Author:-----
% by David L. Jones, Dec-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (nargin < 2), verb = 0; end % default don't show merged PDF in Preview

if isempty(strmatch('.tex',fname(end-3:end)))
   error('FNAME must be a ''*.tex'' LaTeX file!')
end
% ----------------------

% Run latex:
cmd = ['/usr/texbin/pdflatex --shell-escape ' fname];
system(cmd);

% Optionally display typeset PDF:

if (verb>0)
   pdfName = [fname(1:end-4) '.pdf'];
   fprintf('Press any key after LaTeX if finished... \n')
   pause
   eval(['!open -a Preview ' pdfName]);
end
