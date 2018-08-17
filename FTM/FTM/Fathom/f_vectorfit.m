function vects = f_vectorfit(crds,env,scale,offset,labels);
% - plot environmental ordination vectors via multiple linear regression
%
% USAGE: vects = f_vectorfit(crds,env,scale,offset,labels);
%
% crds   = matrix of site scores (rows) x ordination dimensions (cols)
%          e.g., the scores from a PCoA Ordination
%
% env    = matrix of sites (rows) x environmental variables (cols)
%
% scale  = scaling factor for species vectors (default = 1)
%
% offset = label offset                       (default = 0);
%
% labels = cell array of env labels (if empty, autocreate)
%          e.g., labels = {'sal' 'temp' 'depth'};
%
% vects = matrix of site coefficients (rows) x env variables (cols);
%         coefficients define the coordinates of the heads of vectors;
%         vectors are scaled to unit length in FULL dimensional
%         space, so vector length is proportional to correlation b/n each
%         variable & the PLOTTED ordination space.
%
% SEE ALSO: f_biplotEnv2, f_biplotEnv3, f_biplotSpecies, f_pca, f_pcoa, f_nmds

% -----References:-----
% Belbin, L. 1988. PATN Reference Manual. CSIRO Division of Wildlife & Ecology, Canberra.
% Dargie, T.C.D. 1984. On the integrated interpretation of indirect site ordinations:
%   a case study using semi-arid vegetation in southeastern Spain. Vegetatio 55: 37-55.
% Faith, D.P. & R.H. Norris. 1989. Relation of environmental variables with patterns
%   of distribution & abundance of common & rare freshwater macroinvertebrates. Biol.
%   Conserv. 50: 77-98.
% Minchin, P. 1991. DECODA version 2.04. Preliminary Documentation.
%   Australian National University.

% -----Credits:-----
% inspired by vectorfit for "R" by Jari Oksanen<jarioksa@pc112145.oulu.fi>
% http://cc.oulu.fi/~jarioksa/softhelp/vegan.html

% -----Author:-----
% by David L. Jones, June-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 19-July-2001 debugged labels

if (nargin < 3), scale   = 1; end; % set default scaling to 1
if (nargin < 4), offset  = 0; end; % default is no label offset
if (nargin < 5), labels = num2cell([1:size(env,1)]); end; % default vector labels

noCol = size(env,2); % get # of env variables

[yfit,r2,coefs,resid,p] = f_mregress(crds,env,1000); % multiple linear regression

rowsCoefs = size(coefs,1);    % get # of rows in coefs
coefs = coefs(2:rowsCoefs,:); % strip off y-intercept

for i=1:noCol % normalize vectors for each variable to unit length
   vects(:,i) = coefs(:,i)./(sqrt(sum(coefs(:,i).^2))); % make sum-of-squares = 1
end;

vects = vects(1:2,:);
vects = vects';

vects = vects*scale; % scaled vects

hold on;
% ----- overlay vects on PCoA ordination: -----
for j = 1:noCol
   thisVector = [0,0;vects(j,1),vects(j,2)];
   plot(thisVector(:,1),thisVector(:,2),'-k');
   
   if (vects(j,1)>=0 & vects(j,2)>=0)
      deltaX = offset;
      deltaY = offset;
   elseif (vects(j,1)<0 & vects(j,2)>=0)     
      deltaX = -1*offset;
      deltaY =    offset;
   elseif (vects(j,1)<0 & vects(j,2)<0)   
      deltaX = -1*offset;
      deltaY = -1*offset;
   elseif (vects(j,1)>=0 & vects(j,2)<0)
      deltaX =    offset;
      deltaY = -1*offset;
   end;
   
   h = text(thisVector(2,1)+(deltaX*0.5),thisVector(2,2)+(deltaY*0.5),labels(j),'FontSize', 8); % label env vectors
   set(h,'HorizontalAlignment','center');
end;

axis equal;
hold off;

