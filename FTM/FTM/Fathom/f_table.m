function f_table(X,fmt,sep)
% - display a matrix X in the command window
% 
% USAGE: = f_table(X,'fmt','sep')
% 
% X   = input matrix
% fmt = format (default = '%13.6e')
% sep = use spaces ('s'), commas ('c'), or tabs ('t'), to delimit columns
%       (default = 's')
% 
% SEE ALSO: f_copy

% -----References:-----
% After asTable.m by Stefan Baunack<s.baunack@ifw-dresden.de.>, 05-08-99

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), fmt = '%13.6e'; end % default format
if (nargin < 3), sep = 't';     end % default tab delimited

if sep=='s'
   sep = ' ';
elseif sep=='c'
   sep = ',';
elseif sep=='t'
   sep = '\t';
else
   error('SEP must be ''s'', ''c'', or ''t''!')
end
% -------------------------------------

if ~isreal(X), error('X must contain only real elements.'),end
if  ischar(X), error('X is a string matrix.'),end

[nr,nc] = size(X);

% Add leading line break:
disp(' ')

if (nc==1) % row vector:
   % Create output string line by line at once:
   txt = sprintf([fmt '\n'],X(:));

   % Replace +'s with spaces:
   if (sum(fmt=='+')), txt = regexprep(txt,'\+',' '); end 
   
   % Display text:
   disp(txt)

else
   for L = 1:nr
      % Create output string for row L without trailing tab:
      txt = [sprintf([fmt sep],X(L,1:nc-1)) sprintf(fmt,X(L,nc))];

      % Replace +'s with spaces:
      if (sum(fmt=='+')), txt = regexprep(txt,'\+',' '); end
      
      % Display text:
      disp(txt)
   end
end

% Add trailing line break:
disp(' ')

