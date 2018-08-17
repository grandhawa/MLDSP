function result = f_ibc(con,Y,fname,yLabels,sLabels,tol)
% - indicator-species based cluster analysis
% 
% USAGE: result = f_ibc(mst.con,Y,'fname',yLabels,sLabels,tol);
% 
% mst.con = pairs of sites connected by a minimum spanning tree (see f_mst)
% Y       = (transformed) species abundance data
% fname   = base name of file to write results to (e.g., fname = 'output')
% yLabels = cell array of Y labels; if empty, autocreate
%           e.g., yLabels = {'spA' 'spB' 'spC'};
% sLabels = cell array of site labels; if empty, autocreate
%           e.g., sLabels = {'s1' 's2' 's3'};
% tol     = tolerance; the maximum difference between an observed scaled
%           abundance value and 1 before the observed value is no longer
%           considered to be 'relatively high'; values must range from 0-1
%           with the default = 1 SD of the scaled data (for each species)
% 
% result = structure of results with the following fields:
%  .Y    = species abundance values scaled to range from 0 to 1
%  .tol  = tolerance value used for each indicator species
%  .con  = list of connections between sites within tolerance
%  .idxS = index to sites indicative of each indicator species
%  .S    = Boolean values indicating sites indicative of species (1=true, 0=false)
%  .B    = symmetric binary connection matrix
%  .uGrp = unique groups of Boolean values
%  .yTxt = corresponding list of unique groups of indicator species
%  .sTxt = corresponding list of sites where each unique group occurs
% 
% SEE ALSO: f_mst, f_plotNeigh

% -----Notes:-----
% This function performs Indicator-species Based Cluster (IBC) analysis.
% Start with a data set consisting of the abundance and composition of species
% among sampling sites. Create a symmetric dissimilarity matrix among sites
% using 'f_dis' and examine the data using a multivariate ordination
% technique (e.g., 'f_pcoa') to determine which taxa should be considered
% 'indicator species'. Use the same dissimilarity matrix to create a
% minimum spanning tree with 'f_mst' that can be overlaid on the
% ordination. Then, for each indicator species, use the IBC algorithm to:
% (1) search the scaled abundance values to identify the sites where the
% indicator species occurs at relatively high abundances; and (2) trim the
% minimum spanning tree (MST) to only include connections between pairs of
% sites where both have relatively high abundances of the indicator
% species. Those sites comprising the trimmed MST are considered indicative
% of the indicator species and comprise the members of a newly identified
% 'cluster' of sites.
% 
% NOTE: Y is the original (transformed) species abundance date used to
% create the symmetric dissimilarity matrix, that was in turn used to
% construct the MST.
% 
% NOTE: the rows of 'yTxt' and 'sTxt' correspond to the rows of 'uGrp' (i.e., the
% unique groups of indicator species).  

% -----Author:-----
% by David L. Jones, Jan-2013
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
% 
% Mar-2013: now supports yLabels, sLabels, and outputs uGrp, yTxt, and
%           sTxt.
% Apr-2013: don't remove groups with no indicator species; exports tabular
%           reports to files specified in fname;

% -----Check inputs and set defaults:-----
if (nargin < 4), yLabels = []; end % default species labels
if (nargin < 5), tol     = []; end % no tolerance provided

% Check CON:
if (size(con,2)~= 2)
   error('MST.CON must be a list of connections created by F_MST!');
end

% Define files names to write results to:
fnameY = [fname '_y.csv'];
fnameS = [fname '_s.csv'];

% Make sure files doesn't alreay exist:
   if (exist(fnameY,'file')==2)
      error([fnameY ' already exists!'])
   end
   
   if (exist(fnameS,'file')==2)
      error([fnameS ' already exists!'])
   end
   % ----------------------------------------

% Get size of input:
[r,c] = size(Y);

% Scale values relative to column maximum:
Y = Y ./ repmat(max(Y),r,1);

% -----Check Labels:-----
% Autocreate yLabels:
if isempty(yLabels)
   yLabels = cellstr(num2str((1:ncY)'))';
end

% If labels are not cell arrays, try forcing them:
if iscell(yLabels)<1, yLabels = num2cell(yLabels); end

% Make sure labels are of compatible size:
yLabels = yLabels(:);
if size(yLabels,1) ~= c;
   error('Size mismatch of yLables and Y!')
end

% Autocreate sLabels:
if isempty(sLabels)
   sLabels = cellstr(num2str((1:r)'))';
end

% If labels are not cell arrays, try forcing them:
if iscell(sLabels)<1, sLabels = num2cell(sLabels); end

% Make sure labels are of compatible size:
sLabels = sLabels(:);
if size(sLabels,1) ~= r;
   error('Size mismatch of sLables and Y!')
end
% ------------------------


% -----Check tolerance:-----
% Empty?
if isempty(tol)
   tol = (1-std(Y));
end

% One value for all columns?
if numel(tol)==1
   tol = repmat(tol,1,c);
end

% Out of bounds?
if (sum(tol<0) || sum(tol>1)) > 0
   error('TOL must be between 0-1');
end

% Size mismatch?
if numel(tol) ~= size(Y,2)
   error('Size mismatch between Y and TOL!')
end

% Force row vector:
tol = tol(:)';
% --------------------------

% Boolean indicating sites within tolerance:
tol_B = (Y - repmat(tol,r,1)) >= 0; % 1 = within, 0 = below

% Index to connections between sites where BOTH are within tolerance:
S = zeros(r,c); % preallocate
for i = 1:c % each species separately
   % Index to rows of tol_B where sites are within tolerance:
   idxT{i} = find( tol_B(:,i) == 1);
   
   % Index to rows of CON where connections only involve sites within tolerance:
   idxC{i} = find( (ismember(con(:,1),idxT{i}) .* ismember(con(:,2),idxT{i})) == 1);

   % Update list of connections (trimmed version):
   conT{i} = con(idxC{i},:);

   % Index to sites indicative of this species:
   idxS{i} = f_unique(conT{i}(:));
   
   % Boolean indicating sites indicative of this species:
   S(:,i) = ismember([1:r]',idxS{i});
   
   % Create a symmetrical binary connection matrix:
   B{i}      = zeros(r,r);                               % initialize
   idx       = sub2ind([r r],conT{i}(:,1),conT{i}(:,2)); % convert subscript to indices
   B{i}(idx) = 1;                                        % indicate pairs that are connected
   B{i}      = B{i} + B{i}';                             % fill upper tridiagonal
end


% -----Create list of unique groups of species:-----
% Find unique groups of indicator species:
uGrp = f_unique(S);

% Remove groups with no indicator species:
% uGrp(find(sum(uGrp,2)==0),:) = [];

% Generate a table:
yTxt = []; % initialize
for i = 1:size(uGrp,1) % repeat for each unique group
      
   % Make a copy of original labels:
   temp = yLabels;
   
   % Get index to species NOT comprising this group:
   idx = find( uGrp(i,:)==0 );
   
   % Replace missing species with '-'
   temp(idx) = repmat({'-'},1,numel(idx));
   
   % Add indicators species to the list:
   yTxt = [yTxt sprintf('%s, ', temp{:} )];
   yTxt = [yTxt sprintf('\n')];
end

% -----List sites associated with unique groups:-----
% Generate a table:
sTxt = []; % initialize

for i = 1:size(uGrp,1) % repeat for each unique group
   % Get index to sites matching this group:
   idx = find(ismember(S,uGrp(i,:),'rows') == 1);
   txt = sLabels(idx);
   
   % Add sites to the list:
   sTxt = [sTxt sprintf('%s, ', txt{:} )];
   sTxt = [sTxt sprintf('\n')];
end
% -------------------------------------------------------

% Write file:
f_writeTxt(yTxt,[fname '_y.csv']); % list of unique groups of indicator species
f_writeTxt(sTxt,[fname '_s.csv']); % list of sites where each group occurs

% Wrap results up into a structure:
result.Y    = Y;    % scaled abundance values (0-1)
result.tol  = tol;  % tolerance
result.con  = conT; % list of connections between sites within tolerance
result.idxS = idxS; % index to sites indicative of each indicator species
result.S    = S;    % Boolean values indicating sites indicative of species (1=true, 0=false)
result.B    = B;    % symmetric binary connection matrix
result.uGrp = uGrp; % unique groups of Boolean values
result.yTxt = yTxt; % corresponding list of unique groups of indicator species
result.sTxt = sTxt; % list of sites each unique group occurs
