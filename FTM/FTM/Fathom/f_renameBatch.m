function f_renameBatch(cFile,ptn,rep,tag,verb)
% - rename variables across multiple files
%
% USAGE: f_renameBatch(cFile,ptn,rep,tag,verb);
%
% cFile = cell array of names of files to process
% ptn   = regular expression pattern of variable names to match
% rep   = regular expression pattern of variable names to replace
% tag   = tag each file with a variable indicating its name        (default = 1)
% verb  = verbose output of progress                               (default = 1)
%
%  SEE ALSO: f_rename, regexp, saver

% -----Notes:-----
% This function is used to process one or more Matlab '*.mat' files by
% renaming the variables according to a match/replace pattern specified by
% a regular expression (regex). The original file is replaced by one
% containing variables using the new names. However, the original data is
% retained in a backup copy, which has a new filename tagged with the
% current date/time.
%
% This function optionally adds a variable called 'fname', if it doesn't
% previously exist, which specifies the base name of the file and can be
% used by the 'saver' function.

% -----Author:-----
% by David L. Jones, Aug-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin <  4), tag  = 1; end % default tag each file with its name
if (nargin <  5), verb = 1; end % default verbose output of progress

% Check that CFILE is a cell array:
if ~iscell(cFile), error('CFILE must be a cell array'); end

% Check that files in CFILE exist and contain the proper variables:
nFiles = numel(cFile); % get # files

for i = 1:nFiles
   % Check that file exists:
   if (exist(cFile{i},'file')~=2)
      error(['Cannot find file ' cFile{i} '!'])
   end
   
   % Check that specified file is a *.mat file:
   if (~isequal(cFile{1}(end-3:end),'.mat'))
      error(['File ' num2str(i) ' must specify a *.mat file!']);
   end
   
   % Check for presence of specified variables:
   if isempty(who('-file',cFile{i},'-regexp',ptn));
      error(['File ' cFile{i} ' contains no variables matching ' '''' ptn '''' '!']);
   end
   
   % Create time stamp:
   timeStamp = f_clock;
   
   % Create list of files to backup original data to:
   cFileBAK{i} = regexprep(cFile{i},'.mat$',['_' timeStamp '.mat']);
end
% -------------------------------------

% Process each file separately:
for i = 1:nFiles
   
   % Show file being processed:
   if (verb>0)
      if (i==1), fprintf('\n'); end
      fprintf('Processing file: %s...\n',cFile{i});
   end
   
   % Make backup of original file:
   [flag,msg] = copyfile(cFile{i},cFileBAK{i});
   if ~isempty(msg)
      disp();
      disp(msg);
      disp();
      error(['Creating backup of ' cFile{i} ' was NOT successful!'])
   end
   
   % Load old variables into a structure
   S = load(cFile{i});
   
   % Get list of old variable names:
   S_var = fieldnames(S);
   
   % Create list of new variable names:
   S_rep = regexprep(S_var,ptn,rep);
   
   % Rename variables:
   nVar = numel(S_var);
   for j = 1:nVar
      % Generate command string:
      cmd = ['f_renameField(S,' '''' S_var{j} '''' ',' '''' S_rep{j} '''' ');'];
      
      % Evaluate command string:
      eval(cmd);
      
      % Optionally tag each file with a variable indicating its base name:
      if (tag>0) && (~isfield(S,'fname')==1)
         fname   = regexp(cFile{i},'(?<=/)\w+(?=.mat$)', 'match'); % extract base name
         fname   = char(fname{:});                                 % convert cell to char
         S.fname = fname;                                          % add field
      end
   end
   
   % Save all variables in the workspace to a new file:
   eval(['save ' cFile{i} ' -STRUCT S;']);
   
   % Cleanup for next file:
   clear S;
end

% Show completion of batch process:
if (verb>0)
   fprintf('->...DONE!\n');
end
