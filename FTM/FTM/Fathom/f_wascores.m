function scores = f_wascores(config,spp,plt,labels)
% - weighted-averages scores of species for a site ordination
%
% Usage: scores = f_wascores(config,spp,plt,labels);
%
% config  = ordination configuration from nMDS, PCoA etc.
%           (rows = site scores, col = axes)
% spp   = transformed species abundances
%           (rows = sites, col = abundances)
%
% plt   = optionally create plot (default = 1)
%
% labels = optional cell array of species labels (if empty, autocreate);
%           e.g., labels = {'sp1' 'sp2' 'sp3'};
%
% scores = coordinates of species in site-ordination space
%
% SEE ALSO: f_rdaPlot, f_nmds, f_pcoa, f_pca

% -----Notes:-----
% Computing weighted-averages scores for species allows you to derive
% coordinates for species and plot them in the ordination space defined by
% sites. The score for a species along an ordination axis is the average of 
% scores of the sites it occurs in, but weighted by its abundance in each site.
% This allows simultaneous plotting of species and sites in the same ordination
% space, as in RA, DCA, & CCA, etc.

% -----References:-----
% Minchin, P. 1991. DECODA version 2.04. Preliminary Documentation. Australian
%   National University. 
% McCune, B. & M. J. Mefford. 1999. PC-ORD. Multivariate Analysis of Ecological
%   Data, Version 4. MjM Software Design, Gleneden Beach, Oregon, USA. 237 pp.

% -----Credits:-----
% Inspired by "wascores" for "R" by Jari Oksanen<jarioksa@pc112145.oulu.fi>
% http://cc.oulu.fi/~jarioksa/softhelp/vegan.html

% -----Author:-----
% by David L. Jones, July-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Apr-2003: minor changes to plot options
% Apr-2008: updated documentation; commented out unused variables

% -----Check input and set default values:-----
if (nargin < 3), plt   = 1; end;                          % defaults to plot
if (nargin < 4), labels = num2cell([1:size(spp,2)]); end; % default species labels

% If labels are not cell arrays, try forcing them:
if iscell(labels)<1, labels = num2cell(labels); end;

% Check size of inputs:
if (size(config,1) ~= size(spp,1))
   error('# of sites in CONFIG not equal to SPP!');
end
% ---------------------------------------------

% noSites = size(config,1); % get # of sites
noSpp   = size(spp,2);    % get # of species
noAxes  = size(config,2); % get # of axes

scores(noSpp,noAxes) = 0; % preallocate result matrix 

for i = 1:noAxes
   for j = 1:noSpp
      scores(j,i) = weighted_average(config(:,i),spp(:,j));
   end;
end;

% ----- Plot weighted-averages scores w/ ordination:-----
if (plt>0)
   figure; % open new figure window
   hold on;
   box on;
   plot(config(:,1),config(:,2),'b.');
   % plot(scores(:,1),scores(:,2),'rx'); % plot symbols
   for k = 1:noSpp % plot labels
      h = text(scores(k,1),scores(k,2),labels(k)); % plot labels
      set(h,'HorizontalAlignment','center','Color',[1 0 0],'FontSize',8,...
         'FontWeight','bold');
   end;
   hold off;
end;

%%%%%%%%%% SUBFUNCTION: %%%%%%%%%%%%%%%%%%
function wa = weighted_average(x,w)
sumw = sum(w);
wa   = sum(w.*x) / sumw;
% - after meanwt.m from Richard Strauss:
% http://www.biol.ttu.edu/Faculty/FacPages/Strauss/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
