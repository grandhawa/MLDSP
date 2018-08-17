function [zcm,zse,stn] =  f_depthCM(station,net,zMean,zWidth,zStdv,density)
% - 'Depth of Center of Mass' for plankton survey data
%
% USAGE: [zcm,zse,stn] = f_depthCM(station,net,zMean,zWidth,zStdv,density);
%
% station = horizontal sampling station
% net     = vertical sample taken at STATION
% zMean   = mean towing depth of NET (in meters)
% zWidth  = width of depth strata sampled by NET (max - min towing depth);
% zStdv   = standard deviation of towing depth sampled by NET
% density = standardized densities of planktonic taxa collected in NET at
%           STATION; (multiple taxa my be entered column-wise)
%
% stn = station
% zcm = depth of center-of-mass
% zse = standard error
%
% SEE ALSO: f_depthMean

% -----Notes:-----
% This program is used to calculate the vertical 'Depth of Center of Mass'
% at each horizontal sampling station for planktonic taxa collected using 
% vertically stratified sampling gear, such as MOCNESS.
%
% Mean depths are calculated separately for multiple taxa that are entered
% column-wise in DENSITY.

% -----References:-----
% Roepke, A., W. Nellen, & U. Piatkowski. 1993. A comparative study on the
%  influence of the pycnocline on the vertical distribution of fish larvae
%  and cephalopod paralarvae in three ecologically different areas of the
%  Arabian Sea. Deep-Sea Res. II 40(3): 801-819.

% -----Author:-----
% by David L. Jones, May-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if (sum(size(unique([size(station,1) size(net,1) size(zMean,1) size(zWidth,1)...
         size(zStdv,1) size(density,1)]))) ~=2)
   error('All input variables must have same # of rows !');
end

stn = unique(station);    % unique horizontal stations
ns  = size(stn,1);        % # of horizontal stations
nc  = size(density,2);    % # of taxa

zcm(ns,nc) = 0; % preallocate
zse(ns,nc) = 0; 

for t = 1:nc % do for each taxa separately
   
   for k = 1:ns; % find zcm & zse for each station
      idx = find(station == stn(k));
      
      if (sum(density(idx,t))==0)
         zcm(k,t) = 0;
         zse(k,t) = 0;
      else
         idx(find(density(idx,t)==0)) = []; % remove nets with 0 catch
         
         
         n  = size(net(idx),1); % # of nets for station k with catch > 0;
         
         zcm(k,t) = sum(([density(idx,t) .* zWidth(idx)] ./ [repmat(sum(density(idx,t)...
               .* zWidth(idx)),n,1)]) .* zMean(idx));
         zse(k,t) = sqrt(sum(([density(idx,t) .* zWidth(idx)] ./ [repmat(sum(density(idx,t)...
               .* zWidth(idx)),n,1)]).^2 .* (zStdv(idx)).^2))/sqrt(n);   
      end
      
   end
   
end

