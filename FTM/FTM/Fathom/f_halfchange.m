function crds = f_halfchange(config,diss,thres,nthres);
% - scale nMDS configuration to half-change
%
% USAGE: hc = f_halfchange(config,diss,thres,nthres);
%
% config = nMDS configuration (rows = sites, cols = axes)
% diss   = symmetric dissimilarity matrix used for nMDS
%          (rows = cols = # sites)
% thres  = dissimilarity threshold (default = 0.8)
% nthres = min # of values < thres (default = 10)
%
% crds   = coordinates of half-change configuration
%
% SEE ALSO: f_nmds, f_pcoa, f_pca

% -----References:-----
% after Jari Oksanen's "postMDS" from his "Vegan Package for R"
% <jarioksa@pc112145.oulu.fi> http://cc.oulu.fi/~jarioksa/softhelp/vegan.html
% originally based on an unpublished method by Peter Minchin

% -----Author:-----
% ported from "R" to Matlab
% by David L. Jones, July-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set default values:-----
if (size(config,1) ~= size(diss,1)),
   error('CONFIG & DISS do not have compatible sizes!');
end;

if (f_issymdis(diss) == 0)
   error('Input DISS must be a square symmetric distance matrix');
end;

if (nargin < 3), thres  = 0.8; end; % set default threshold
if (nargin < 4), nthres = 10;  end; % set default threshold
% -------------------------------------------

diss   = f_unwrap(diss); % unwrap lower tri-diagonal as vector
ord    = f_unwrap(f_euclid(config')); % ditto

diss = diss(find(diss<thres)); % extract elements from diss that are < threshold
ord  = ord(find(diss<thres));  % extract corresponding elements from ord

% regress dissimilarities on distances:
[yfit,r2,coefs,resid,p] = f_mregress(diss,ord);

hc = (1-coefs(1))/2/coefs(2);

crds = config/hc;

