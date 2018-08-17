function [res,resLabels] = f_bioenv(yDis,x,labels,rank,metric,trim,out)
% - correlation between distance matrix and all subsets of X (BEST or BIOENV)
%
% USAGE: [res,resLabels] = f_bioenv(yDis,x,labels,rank,metric,trim,out);
%
% yDis    = symmetric distance matrix  
% x       = matrix to subset (rows = obs, cols = variables)
%
% labels = cell array of variable labels x
%          e.g., labels = {'temp' 'sal' 'depth' 'O2'}';
%
% metric = distance metric to construct xDis
%          0 = Euclidean (default); 1 = Bray-Curtis
%
% rank   = rank correlation (default = 1)
%
% trim   = return only this many of the top Rho's per subset size class
%          (0 = return all, default)
%
% out    = send results to screen (= 1, default)
%          or cell array with filename; e.g., out = {'results.txt'}
%          NOTE: existing file with same name will be DELETED!
%
% res       = cell array, 1st col is Rho, 2nd:end are variable indices
% resLabels = cell array of variable names 
% Tabulated results are also sent to screen or file, depending on 'out'
%
% Note: yDis and X must have same # of rows
% 
% SEE ALSO: f_anosim, f_mantel

% -----Notes:-----
% Note the equivalent function in Clarke's PRIMER program, previously known
% as BIOENV, is now called BEST.

% -----References:-----
% Clarke, K. R. 1993. Non-parametric multivariate analyses of changes
% in community structure. Aust. J. Ecol. 18: 117-143.
%
% Clarke,K. R. & M. Ainsworth. 1993. A method of linking multivariate community
% structure to environmental variables. Mar. Ecol. Prog. Ser. 92:205-219.
%
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
% Elsevier Science BV, Amsterdam. xv + 853 pp.

% -----Dependencies:-----
% combnk.m STATISTICS Toolbox (could be replaced by choosenk)

% -----Author(s):-----
% by David L. Jones, Feb-2001
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2003: rank correlation now optional; improved sort accounts for
%           negative correlations, changed input format (colums
%           are variables, row are obserations)
% Mar-2010: replaced '&' with '&&'
% Jan-2011: updated documentation, force labels to col vs. row vector

% -----Set defaults & check input:-----
if (nargin < 4), rank   = 1; end; % rank correlation by default
if (nargin < 5), metric = 0; end; % set default metric to Euclidean
if (nargin < 6), trim   = 0; end; % don't trim results by default
if (nargin < 7), out    = 1; end; % send results to screen by default

if (f_issymdis(yDis) == 0)
   error('yDIS must be square symmetric distance matrix !');
end;

% Check that yDis & x have equal row count:
if size(yDis,1) ~= size(x,1)
   error('YDIS and X must have same # of rows !');
end;

% Force cell array of labels to be col vs. row vector
labels = labels(:);
% -------------------------------------

% If 'labels' are not a cell array, try forcing:
if iscell(labels)<1, labels = num2cell(labels); end;

% Setup destination of tabulated results:
if (iscell(out)<1) % send to screen
   fid = 1;
else % open file for writing:
   fid = fopen(char(out{:}),'w'); % char needed for fopen
end;
fidStr = num2str(fid); % string needed for eval()

% Unwrap lower tridiagonal of yDis:
pri_vect = f_unwrap(yDis);

% Standardize variables of x:
if (metric == 0), x = f_stnd(x); end;

noVar  = size(x,2);   % # of variables in X
noComb = 2^noVar - 1; % # of possible subsets

fprintf('\nThere are %d possible subsets of %d variables \n', noComb, noVar);

for i = 1:noVar %%--do for each size "class" of subsets--%%
   subsets     = combnk([1:noVar],i); % all possible subsets of size i
   noSubsets   = size(subsets,1);     % # of subsets (rows) of size i
   rho(noSubsets) = 0;                % preallocate results array
   
   fprintf('  Processing %4d subsets of %3d variables \n', noSubsets, i);
   
   for j = 1:noSubsets %%--do for each subset of size i--%%
      sMatrix = x(:,[subsets(j,:)]); % extract all rows of x, but only these cols
      
      if (metric==0)
         sec_dist = f_euclid(sMatrix');     % Euclidean distance
      else
         sec_dist = f_braycurtis(sMatrix'); % Bray-Curtis metric
      end;
      
      sec_vect = f_unwrap(sec_dist);             % unwrap lower tridiagonal  
      rho(j)   = f_corr(pri_vect,sec_vect,rank); % matrix correlation
   end;
   
   % Sort results by descending abs(rho):
   res{i} = flipud(sortrows([abs(rho') rho' subsets],1));
   res{i}(:,1) = []; % trim abs(rho')
   
   % Optionally trim results:
   if (trim>0) && ((i>1) && (size(res{i},1)>trim)), res{i} = res{i}(1:trim,:); end;   
   
   % Extract variable labels:
   resLabels{i} = labels([res{i}(:,2:end)]);
   
   % Cleanup before next iteration:
   rho = [];
end;

% Transpose last iterations (col vector to row vector);
resLabels{noVar} = resLabels{noVar}';


%-----Output results to screen or file:-----
fprintf(fid,'\n ========================================== \n');
fprintf(fid,' Rho    Variables \n');
fprintf(fid,' ========================================== \n');
labelStr = ''; % initialize variable
for k = 1:noVar %%--do for each subset size class--%%
   labelStr = [labelStr ' %s'];
   noRho = size(res{k},1); % # of Rho's saved in res cell array
   eval(['fprintf(' fidStr ',' '''\n %d\n'',' num2str(k) ');']);
   for m = 1:noRho 
      eval(['fprintf(' fidStr ',' '''%7.4f' labelStr '\n'', res{k}(m,1),resLabels{k}{m,:});']);
   end;
end;
fprintf(fid,'\n'); % terminating linefeed
if fid ~= 1
   status = fclose(fid); % close output file for writing
   fprintf('\n Done!...Results saved to file: %s \n', char(out{:}));
end;

