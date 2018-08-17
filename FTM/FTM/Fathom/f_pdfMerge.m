function f_pdfMerge(fnames,pdfName,crop,verb)
% - merge a series of PS (or PDF) files into a single PDF
%
% USAGE: f_pdfMerge(fnames,'pdfName',crop,verb)
%
% fnames  = cell array specifying files to merge; e.g., {'plot_1.pdf' 'plot_2.pdf'};
%           - use '*.ps' or '*.pdf' to merge all in the current directory
% pdfName = name of output file (e.g., = 'merged.pdf')
% crop    = crop PDF to bounding box                               (default = 1)
% verb    = display merged PDF in Preview.app                      (default = 1)
%
%
% SEE ALSO: f_pdf, f_pdfCrop, f_pdfDistill, f_pdfLatex, f_pdfSplit

% -----Author:-----
% by David L. Jones, Nov-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% ----- Notes: -----
% This function has been developed under 'MacOS X' and requires Ghostscript and
% pdfcrop, the latter is a Perl script available from http://www.ctan.org and
% should be placed in a directory that is in your user path (e.g., ~/bin).
%
% CROP = 1: this is especially useful when merging PS files, since Ghostscript
% will create a merged PDF with page sizes equal to the default for our
% installation, which is usually letter or A4. Cropping will remove the white
% border around each page.

% -----Check input & set defaults:-----
if (nargin < 3), crop = 1; end % default crop each merged PDF
if (nargin < 4), verb = 1; end % default show merged PDF in Preview

% Delete output file if it already exits:
x = size(dir(pdfName),1);
if (x>0), delete(pdfName); end;
% -------------------------------------


% Generate list of *.PDF or *.PS files in the current directory
if  ischar(fnames) && strcmp(fnames,'*.ps') || ischar(fnames) && strcmp(fnames,'*.pdf')
   dirVar = dir(fnames);
   for i = 1:size(dirVar,1)
      files(i) = {dirVar(i).name};
   end
else
   files = fnames(:);                    % use files specified on command line
end

% Append all filenames to a single variable:
files    = files(:); % force as columns
n        = size(files,1);
inputVar = ''; % initialize
for i = 1:n
   inputVar = [inputVar ' ' files{i}];
end

% Create Ghostscript command:
cmdVar = ['!/usr/local/bin/gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite '...
   '-sOutputFile=' pdfName];

% Run external OS command:
eval([cmdVar ' ' inputVar]);

% Optionally crop each page to the bounding box:
if (crop>0)
   f_pdfCrop(pdfName,pdfName)
end

% Optionally display merged PDF:
if (verb>0)
   if (crop>0)
      fprintf('\nPress any key after the PDF has been cropped... \n')
      pause
   end
   
   eval(['!open -a Preview ' pdfName]);
end
