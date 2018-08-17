function f_pdf(fname,noPDF,verb,level)
% - export current figure to PDF file
%
% USAGE: f_pdf('fname',noPDF,verb,level)
%
% fname = filename to export to, WITHOUT the extension.
% noPDF = generate PostScript only, no PDF              (default = 1)
% verb  = display PDF in Preview.app                    (default = 1)
% level = use Postscript Level 1 vs. 2                  (default = 1)
%
% SEE ALSO: f_pdfMerge, f_pdfDistill

% ----- Notes: -----
% This function give similar results as PRINT -DPDF, but produces a PDF file
% that is cropped to the figure size. This is useful if you're creating PDF
% graphics for inclusion in page-layout programs, etc.
%
% Acrobat Distiller will use the default SETTINGS when creating the PDF.
%
% This function was developed under 'Mac OS X'
% 

% -----Supported Postscript Fonts:-----
% The Matlab Postscript/Ghostscript drivers only support the following fonts; if
% you select anything else, it will be replaced by 'Courier':
% 
% 'AvantGarde', 'Helvetica-Narrow', 'Times-Roman', 'Bookman', 'NewCenturySchlbk'
%  'ZapfChancery', 'Courier', 'Palatino', 'ZapfDingbats', 'Helvetica', 'Symbol'

% -----Author:-----
% by David L. Jones, Dec-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2003: improved documentation, added noPDF option
% Apr-2003: bug fix
% Aug-2006: changed to PS level 1 (from 2), works better
% Jun-2007: edited to work under 'Mac OS X'
% Nov-2008: added support for 'verb' and Acrobat Distiller
% Jan-2009: added option for PS Level 2
% Feb-2009: added documentation
% Oct-2009: corrected 'nargin' for 'level'
% Jun-2012: now check operating system
% Oct-2012: changed psName from '.ps' to '.eps' so graphic will display in
%           Mac OS X's by tapping the spacebar

% -----Check input & set defaults:-----
if (nargin < 2), noPDF = 1; end % don't generate PDF by default
if (nargin < 3), verb  = 1; end % display output by default
if (nargin < 4), level = 1; end % use postscript level 1 by default

if (ismac<1)
   error('This function requires MacOS X!')
end
% -------------------------------------

psName  = [fname '.eps'];
pdfName = [fname '.pdf'];

if (level==1)
   % EPS level 1 color:
   eval(['print -depsc ' psName]);
else
   % EPS level 2 color:
   eval(['print -depsc2 ' psName]);
end

% Optionally create a PDF:
if (noPDF<1)
   % Create PDF with Distiller:
   eval(['!open -a ''Acrobat Distiller''' ' ' psName]);
   fprintf('\nPress any key when Distiller is finished... \n');
   pause
   delete(psName)
end

% Optionally display PS or PDF in Preview.app:
if (noPDF<1) && (verb>0)
   eval(['!open -a Preview ' pdfName]);
elseif (noPDF==1) && (verb>0)
   eval(['!open -a Preview ' psName]);
end
