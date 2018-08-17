function pdfNames = f_pdfSplit(fname)
% - split a multipage PDF file into separate single-page PDF's
% 
% USAGE = f_pdfSplit('fname')
% 
% fname    = name of input file (e.g., = 'merged.pdf')
% pdfNames = cell array listing names of extracted PDF files
% 
% SEE ALSO: f_pdf, f_pdfCrop, f_pdfDistill, f_pdfLatex, f_pdfMerge


% -----Author:-----
% by David L. Jones, Dec-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Use Ghostscript to query the number of pages in PDF:
cmd = ['/usr/local/bin/gs -dNODISPLAY -sPDFname=' fname ' -sOutputFile=- '...
   '~/bin/pdfpagecount.ps'];
[null,result] = system(cmd);

% Get index corresponding to last digit in GS output:
idx = regexp(result,'\d');

% Convert from string to number:
n = str2num(result(idx(end)));

% Call GS to extract each page as a separate file:
for i = 1:n
   % Generate list of output files as a cell array:
   pdfNames{i} = [fname(1:end-4) '_' num2str(i,'%03d') '.pdf'];
   
   % Call GS:
   cmd = ['/usr/local/bin/gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH '...
      '-dFirstPage=' num2str(i) ' -dLastPage=' num2str(i) ' -sOutputFile='...
      pdfNames{i} ' ' fname];
system(cmd);
end

% Make row vector:
pdfNames = pdfNames';





