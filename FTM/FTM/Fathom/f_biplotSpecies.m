function [biplot,Rsq] = f_biplotSpecies(crds,species,rank,iter,scale,offset,sLabels,wt)
% - create species vectors for ordination distance biplot
%
% Usage: [biplot,rsq] = f_biplotSpecies(crds,species,rank,iter,scale,offset,sLabels,wt);
%
% ----- Input: -----
% crds  = matrix of ordination coordinates
%           (rows = sites; cols = eigenvectors or dimensions)
%
% species = matrix of (transformed) species abundances
%           (rows = sites, cols = species)
%
% rank = type of correlation (0 = Pearson's [default], 1 = Spearman's)
%
% iter    = number of iterations for randomized probabilities (default = 0)
%
% scale   = scaling factor for species vectors (default = 1)
%
% offset  = label offset (default = 0);
%
% sLabels = cell array of species labels (if empty, autocreate)
%           e.g., sLabels = {'sp1' 'sp2' 'sp3'};
%
% wt     = weight species scores (default = 1)
%
% ----- Output: -----
% biplot  = 2 column matrix specifying x-y coordinates of endpoints
%           for {scaled} species vectors to overlay on ordination
%
% Rsq     = column matrix of correlation with each axis
%           (randomized significance provided when iter>0)
%
% SEE ALSO: f_biplotEnv2, f_biplotEnv3, f_vectorfit

% ----- References: --------------------------------------------------
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp.
% Legendre, P. and E. Gallagher. Ecologically meaningful transformations
%   for ordination biplots of species data. Oecology 129: 271-280.
% ---------------------------------------------------------------------

% -----Author:-----
% by David L. Jones, April-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 30-Mar-02: calls f_corr vs. corr
% Apr-2003: added wt as an option
% Oct-2010: replaced & with &&

% ----- Check input and set default values: -----
if (size(crds,1) == size(species,1))<1
   error('The # of rows (sites) in SPECIES & crds must be equal!');
end;

if (nargin < 3), rank    = 0; end; % set default correlation to Pearson's
if (nargin < 4), iter    = 0; end; % set default iterations to 0
if (nargin < 5), scale   = 1; end; % set default scaling to 1
if (nargin < 6), offset  = 0; end; % default is no label offset
if (nargin < 7), sLabels = num2cell([1:size(species,1)]); end; % default species labels
if (nargin < 8), wt      = 1; end; % weight species scores by default


% if labels are not cell arrays, try forcing them:
if iscell(sLabels)<1, sLabels = num2cell(sLabels); end;
% -----------------------------------------------

% get 1st 2 axes for 2-d biplot:
x = crds(:,1); y = crds(:,2);

nloop = size(species,2); % get # of species
biplot(nloop,2) = 0;     % preallocate 2-d array


% get standard deviations:
x_std = std(x);
y_std = std(y);

for i = 1:nloop
   thisSpecies = species(:,i); % extract each species separately
   species_std = std(thisSpecies);
   
   [xR,xP] = f_corr(thisSpecies,x,rank,iter);  % correlation with 1st axis
   [yR,yP] = f_corr(thisSpecies,y,rank,iter);  % correlation with 2nd axis
   
   if (wt>0)
      biplot(i,1) = xR*species_std/x_std;
      biplot(i,2) = yR*species_std/y_std;
   else
      biplot(i,1) = xR;
      biplot(i,2) = yR;
   end
   
   Rsq(i,1)    = xR.^2*100; % actual R-squared correlation
   Rsq(i,2)    = yR.^2*100;
   if iter>0
      Rsq(i,3)    = xP; % randomized probability
      Rsq(i,4)    = yP;
   end;
end;   

biplot = biplot*scale; % scaled biplot

hold on;
% ----- plot the vectors: -----
for j = 1:nloop
   thisVector = [0,0;biplot(j,1),biplot(j,2)];
   plot(thisVector(:,1),thisVector(:,2),'-k');
   
   if (biplot(j,1)>=0 && biplot(j,2)>=0)
      deltaX = offset;
      deltaY = offset;
   elseif (biplot(j,1)<0 && biplot(j,2)>=0)     
      deltaX = -1*offset;
      deltaY =    offset;
   elseif (biplot(j,1)<0 && biplot(j,2)<0)   
      deltaX = -1*offset;
      deltaY = -1*offset;
   elseif (biplot(j,1)>=0 && biplot(j,2)<0)
      deltaX =    offset;
      deltaY = -1*offset;
   end;
   
   h = text(thisVector(2,1)+(deltaX*0.5),thisVector(2,2)+(deltaY*0.5),sLabels(j)); % label species vectors
   set(h,'HorizontalAlignment','center');
end;

axis equal;
hold off;