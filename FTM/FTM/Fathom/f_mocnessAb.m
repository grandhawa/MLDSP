function A = f_mocnessAb(Y,DI,V,tow)
% - returns abundance (# under 10 m^2 of sea surface)
%
% USAGE: A = f_mocnessAb(n,DI,V,trawl);
% Y   = column vector of absolute counts of a single taxon in each net
% DI  = depth interval (m)
% V   = volume filtered for each net
% tow = column vector of integers specifying a unique MOCNESS tow

%  -----Author:-----
% by David L. Jones,, Oct-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

uGrp = f_unique(tow); % unique groups, unsorted
nGrp = size(uGrp,1);  % # of groups 
A    = nan(nGrp,1);   % preallocate

% Make intermediate calculation for each net:
YDIV = Y .* (DI./V);

% Group sum within each tow:
for i=1:nGrp
   idx = find(tow == uGRP(i)); % get indices for members group i
   A(i)= sum(YDIV(idx));       % group sum
end
