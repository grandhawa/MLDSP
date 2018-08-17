function [md_all,sd_all,md_stn,stn] = f_depthMean(station,net,zMean,density)
% - weighted 'Mean Depth of Occurrence' for plankton survey data
%
% USAGE: [md_all,sd_all,md_stn,stn] = f_depthMean(station,net,zMean,density);
%
% station = horizontal sampling stations
% net     = vertical sample taken at STATION
% zMean   = mean towing depth of NET (in meters)
% density = standardized densities of planktonic taxa collected in NET at
%           STATION; (multiple taxa may be entered column-wise)
%
% md_all = mean depth of occurrence (all stations combined)
% sd_all = standard deviation       (all stations combined)
% md_stn = mean depth of occurrence (each station separately)
% stn    = station
%
% SEE ALSO: f_depthCM

% -----Notes:-----
% This function is used to calculate a weighted 'Mean Depth of Occurrence'
% for planktonic taxa collected using vertically stratified sampling gear,
% such as MOCNESS. Mean Depth's are calculated separately for all stations
% combined (MD_ALL) and for each horzontal sampling station (MD_STN).
%
% Mean depths are calculated separately for multiple taxa that are entered
% column-wise in DENSITY.

% -----References:-----
% Brodeur, R. D., and W. C. Rugen. 1994. Diel vertical distribution of
% ichthyoplankton in the northern Gulf of Alaska. Fish. Bull. (US) 92(2):
% 223-235.

% -----Author:-----
% by David L. Jones, May-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (sum(size(unique([size(station,1) size(net,1) size(zMean,1) size(density,1)]))) ~=2)
   error('All input variables must have same # of rows !');
end

stn = unique(station); % unique horizontal stations
ns  = size(stn,1);     % # of horizontal stations
nc  = size(density,2); % # of taxa

md_all(1,nc)      = 0; % preallocate
sd_all(1,nc)      = 0;
md_stn(ns,nc)     = 0;
density_stn(ns,1) = 0;


for t = 1:nc % do for each taxon separately
   
   % -----Weighted Mean Depth of Occurrence:-----
   
   % Mean Depth for each station separately:
   for k = 1:ns; 
      idx = find(station == stn(k)); % extract nets for each stn separately
            
      if (sum(density(idx,t))==0)
         md_stn(k,t) = 0; % no catch at stn k for taxon t
      else
         numer = sum(density(idx,t) .* zMean(idx));
         denom = sum(density(idx,t));
         md_stn(k,t) = sum(numer)/sum(denom); % Mean depth at stn k for taxon t
      end
   end
   
   
   % Mean Depth for all stations:
   for k = 1:ns; 
      idx = find(station == stn(k)); % extract nets for each stn separately
      
      numer(k) = sum(density(idx,t) .* zMean(idx));
      denom(k) = sum(density(idx,t));
   end
   
   md_all(t) = sum(numer)/sum(denom); % Mean depth at all stations for taxon t
   
   
   % -----Standard Deviation:-----
   
   % Total density of larvae at each station:
   for k = 1:ns; 
      idx = find(station == stn(k));        % extract nets for each stn separately
      density_stn(k) = sum(density(idx,t)); % total density at stn k for taxon t
   end
   
   sd_all(t) = sqrt((ns/[(sum(density_stn)^2) * (ns-1)]) *  (sum([(density_stn).^2] .*...
      [(md_stn(:,t) - md_all(t)).^2])));
   
end




