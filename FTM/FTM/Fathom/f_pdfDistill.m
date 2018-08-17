function pdfNames = f_pdfDistill(psNames)
% - distill Postscript files to PDF using Acrobat Distiller
%
% USAGE: f_pdfDistill(psNames)
%
% psNames  = cell array listing files to distill; e.g., {'plot_1.ps' 'plot_2.ps'};
%           - use 'dir' to distill all *.ps files in the current directory
%
% pdfNames = cell array listing names of distilled PDF files
%
% SEE ALSO: f_pdf, f_pdfCrop, f_pdfLatex, f_pdfMerge, f_pdfSplit

% ----- Notes: -----
% This function has been developed under 'MacOS X' and requires Adobe Acrobat
% Distiller. Ghostscript could be used but it doesn't create PDF's with page
% sizes cropped to the size of the bounding box of the image. Rather pages are
% sized according to the default of your Ghostscript installation, which usually
% results in a white page border placed around each image.
% 
% - Note that you can distill a PS file using Ghostscript using the 'f_pdfMerge'
% function, specifying only 1 file as input. Also, 'f_pdf' uses Apple's
% Preview.app to distill a PS file created from the current figure window.

% -----Author:-----
% by David L. Jones, Nov-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2009: improved documentation

% Generate list of input files as a cell array:
if (ischar(psNames) && strcmp(psNames,'dir')) % list PS's in the current directory
   dirVar = dir('*.ps');
   for i = 1:size(dirVar,1)
      files(i) = {dirVar(i).name};
   end
else
   files = psNames(:); % use files specified on command line
end

% Append all filenames to a single variable:
files    = files(:); % force as columns
n        = size(files,1);
inputVar = ''; % initialize
for i = 1:n
   inputVar = [inputVar ' ' files{i}];
end

fprintf('\nPlease wait... \n');
for i = 1:n
   % Distill to PDF:
   eval(['!open -a ''Acrobat Distiller''' ' ' files{i}]);
end
fprintf('\nPress any key when Acrobat Distiller is finished... \n');
pause

% Replace '*.ps' with *.pdf':
pdfNames = regexprep(files,'.ps','.pdf');


