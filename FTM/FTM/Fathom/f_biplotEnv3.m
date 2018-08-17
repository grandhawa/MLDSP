function [biplot,Rsq] = f_biplotEnv3(crds,env,special,iter,scale,offset,sLabels,plotflag,minP);
% - create environmental vectors for 3-d nMDS ordination distance biplot
%
% Usage: [biplot,Rsq] = f_biplotEnv3(crds,env,{special},{iter},{scale},{offset},{sLabels},{plotflag},{minP});
%
% ----- Input: -----
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
% plotflag = plot axis 1 vs. 2 (0 = default) or 1 vs. 3 (1)
%
% minP = if p-value of correlation is > minP, then correlation is NOT used (default = 0.05); 
% 
% ----- Output: -----
% biplot  = 2 column matrix sepcifying x-y coordinates of endpoints
%           for {scaled} env vectors to overlay on ordination
%
% Rsq     = column matrix of correlation with each axis
%           (randomized significance provided when iter>0)
%
% See also: f_biplotEnv2, f_biplotSpecies, & f_vectorfit

% -----Author(s):-----
% by David L. Jones, July-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 31-Jul-01: modified for nMDS instead of PCoA
% 30-Mar-02: calls f_corr vs. corr

% ----- References: --------------------------------------------------
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp. [page 586]
% Legendre, P. (personal communication)
% ---------------------------------------------------------------------

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
if (nargin < 8), plotflag = 0; end; % plot axis 1 vs. 2
if (nargin < 9), minP = 0.05; end; % default pLevel

% if labels are not cell arrays, try forcing them:
if iscell(sLabels)<1, sLabels = num2cell(sLabels); end;

% get 3 dimensions of ordination plot for 2-d biplot:
x = crds(:,1); y = crds(:,2); z = crds(:,3);

nloop = size(env,2); % get # of env
biplot(nloop,3) = 0; % preallocate 2-d array


% get scaling factor for each axis:
[scores,evals,expl] = f_pca(crds,0);

x_factor = sqrt(expl(1)); % Legendre & Legendre, 1998 page 586
y_factor = sqrt(expl(2));
z_factor = sqrt(expl(3));

for i = 1:nloop
   thisEnv = env(:,i); % extract each Env separately
   
   [xR,xP] = f_corr(thisEnv,x,special,iter);  % correlation with 1st axis
   if xP > minP, xR = 0; end;
   [yR,yP] = f_corr(thisEnv,y,special,iter);  % correlation with 2nd axis
   if yP > minP, yR = 0; end;
   [zR,zP] = f_corr(thisEnv,z,special,iter);  % correlation with 3rd axis
   if zP > minP, zR = 0; end;
   
   
   % scale according to "% Variation Explained" by each axis:
   biplot(i,1) = xR*x_factor; 
   biplot(i,2) = yR*y_factor;
   biplot(i,3) = zR*z_factor;
   
   
   Rsq(i,1)    = xR.^2*100; % actual R-squared correlation
   Rsq(i,2)    = yR.^2*100;
   Rsq(i,3)    = zR.^2*100;
   
   if iter>0
      Rsq(i,4)    = xP; % randomized probability
      Rsq(i,5)    = yP;
      Rsq(i,6)    = zP;
   end;
end;   

biplot = biplot*scale; % scaled biplot

figure;
hold on;

if (plotflag == 0)   
   % ----- plot the vectors for axis 1-2: -----
   for j = 1:nloop
      thisVector = [0,0;biplot(j,1),biplot(j,2)];
      plot(thisVector(:,1),thisVector(:,2),'-k');
      
      if (biplot(j,1)>=0 & biplot(j,2)>=0)
         deltaX = offset;
         deltaY = offset;
      elseif (biplot(j,1)<0 & biplot(j,2)>=0)     
         deltaX = -1*offset;
         deltaY =    offset;
      elseif (biplot(j,1)<0 & biplot(j,2)<0)   
         deltaX = -1*offset;
         deltaY = -1*offset;
      elseif (biplot(j,1)>=0 & biplot(j,2)<0)
         deltaX =    offset;
         deltaY = -1*offset;
      end;
      
      h = text(thisVector(2,1)+(deltaX*0.5),thisVector(2,2)+(deltaY*0.5),sLabels(j)); % label env vectors
      set(h,'HorizontalAlignment','center');
   end;
   xlabel('Axis 1'); ylabel('Axis 2'); 
else 
   % ----- plot the vectors for axis 1-3: -----
   for j = 1:nloop
      thisVector = [0,0;biplot(j,1),biplot(j,3)];
      plot(thisVector(:,1),thisVector(:,2),'-k');
      
      if (biplot(j,1)>=0 & biplot(j,3)>=0)
         deltaX = offset;
         deltaZ = offset;
      elseif (biplot(j,1)<0 & biplot(j,3)>=0)     
         deltaX = -1*offset;
         deltaZ =    offset;
      elseif (biplot(j,1)<0 & biplot(j,3)<0)   
         deltaX = -1*offset;
         deltaZ = -1*offset;
      elseif (biplot(j,1)>=0 & biplot(j,3)<0)
         deltaX =    offset;
         deltaZ = -1*offset;
      end;
      
      h = text(thisVector(2,1)+(deltaX*0.5),thisVector(2,2)+(deltaZ*0.5),sLabels(j)); % label env vectors
      set(h,'HorizontalAlignment','center');
   end;
   xlabel('Axis 1'); ylabel('Axis 3'); 
end;
axis equal;
hold off;