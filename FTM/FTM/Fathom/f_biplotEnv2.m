function [biplot,Rsq] = f_biplotEnv2(crds,env,special,iter,scale,offset,sLabels)
% - create environmental vectors for 2-d nMDS ordination distance biplot
%
% Usage: [biplot,Rsq] = f_biplotEnv2(crds,env,{special},{iter},{scale},{offset},{sLabels});
%
% crds   = matrix of nMDS coordinates (rows = sites; cols = dimensions)
%
% env    = matrix of (transformed) environmental variables
%          (rows = sites, cols = variables)
%
% special = type of correlation (0 = Pearson's, 1 = Spearman's [default])
%
% iter    = number of iterations for randomized probabilities (default = 0)
%
% scale   = scaling factor for env vectors (default = 1)
%
% offset  = label offset (default = 0);
%
% sLabels = cell array of vector labels (if empty, autocreate)
%           e.g., sLabels = {'sal' 'tmp' 'elev'};
%
% ----- Output: -----
% biplot  = 2 column matrix sepcifying x-y coordinates of endpoints
%           for {scaled} env vectors to overlay on ordination
%
% Rsq     = column matrix of correlation with each axis
%           (randomized significance provided when iter>0)
%
% SEE ALSO: f_biplotEnv3, f_biplotSpecies, & f_vectorfit

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp. [page 586]
% Legendre, P. 2001. (personal communication)

% -----Author(s):-----
% by David L. Jones <djones@marine.usf.edu>, July-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 31-Jul-01: modified for nMDS instead of PCoA
% 30-Mar-02: calls f_corr vs. corr
% Nov-2003: made compatible with f_rdaPlot
% Oct-2010: replaced & with &&


% ----- Check input and set default values: -----
if (size(crds,1) == size(env,1))<1
   error('The # of rows (sites) in ENV & crds must be equal!');
end;

if (nargin < 2)
   error('At least 2 input variables are required!');
end;

if (nargin < 3), special = 1; end; % set default correlation to Pearson's
if (nargin < 4), iter    = 0; end; % set default iterations to 0
if (nargin < 5), scale   = 1; end; % set default scaling to 1
if (nargin < 6), offset  = 0; end; % default is no label offset
if (nargin < 7), sLabels = num2cell([1:size(env,1)]); end;   % default vector labels

% If labels are not cell arrays, try forcing them:
if iscell(sLabels)<1, sLabels = num2cell(sLabels); end;

% Get 1st 2 dimensions of ordination plot for 2-d biplot:
x = crds(:,1); y = crds(:,2);

nloop = size(env,2); % get # of env
biplot(nloop,2) = 0; % preallocate 2-d array

% Get scaling factor for each axis:
[scores,evals,expl] = f_pca(crds,0,1);

x_factor = sqrt(expl(1)); % Legendre & Legendre, 1998 page 586
y_factor = sqrt(expl(2));

for i = 1:nloop
   thisEnv = env(:,i); % extract each Env separately
   
   [xR,xP] = f_corr(thisEnv,x,special,iter);  % correlation with 1st PCoA axis
   [yR,yP] = f_corr(thisEnv,y,special,iter);  % correlation with 2nd PCoA axis
   
   % Scale according to "% Variation Explained" by each axis:
   biplot(i,1) = xR*x_factor; 
   biplot(i,2) = yR*y_factor;
   
   Rsq(i,1)    = xR.^2*100; % actual R-squared correlation
   Rsq(i,2)    = yR.^2*100;
   if iter>0
      Rsq(i,3)    = xP; % randomized probability
      Rsq(i,4)    = yP;
   end
end   

biplot = biplot*scale; % scaled biplot

hold on;
% ----- plot the vectors: -----
for j = 1:nloop
   %    thisVector = [0,0;biplot(j,1),biplot(j,2)];
   %    plot(thisVector(:,1),thisVector(:,2),'-k');
   
   f_arrow([0 0],[biplot(j,1) biplot(j,2)],'size',[0.125]*0.75,'angle',20,'Color','b')
   
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
   
   h = text(biplot(j,1)+(deltaX*0.5),biplot(j,2)+(deltaY*0.5),sLabels(j)); % label env vectors
   set(h,'FontSize',8,'HorizontalAlignment','center','Color','b');
end;

axis equal;
hold off;