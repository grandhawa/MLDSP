function result = f_NSIDCimport(fname)
% - import a National Snow & Ice Data Center (NSIDC) sea ice coverage binary file
% 
% USAGE: result = f_NSIDCimport('fname',lonlat);
% 
% fname  = file name of NSIDC binary file
% 
% result = structure of results with the following fields:
%   .txt = text header
%   .P   = percent of sea ice coverage (0 - 10%)
% 
% SEE ALSO: f_NSIDCinterp
 
% -----Notes:-----
% National Snow and Ice Data Center (NSIDC)
% 
% http://nsidc.org/data/docs/daac/nsidc0081_ssmi_nrt_seaice.gd.html
% 
% Data format: One-byte flat scaled binary integers (preceded by a 300-byte header);
% a PNG image is provided with each data file. 
% 
% Spatial coverage: North and south polar regions; 25 km resolution 
% 
% Grid type/size: See Polar Stereographic Projections and Grids. Grid size
% varies by region: North: 304 columns x 448 rows South: 316 columns x 332
% rows
% 
% Data range: 0 to 250 (fractional coverage scaled by 250; percentage
% scaled by 2.5) 
% 
% File naming convention: nt_yyyymmdd_fxx_nrt_R.bin
% 
% File size: North: 136492 bytes
%            South: 105212 bytes 
% 
% Parameter: 	Sea ice concentration FTP:
% http://nsidc.org/forms/nsidc-0081_or.html
% 
% Projection: polar stereographic projection
% 
% Parameter Range: The sea ice concentration floating-point values
% (fractional coverage ranging from 0.0 to 1.0) are multiplied by a scaling
% factor of 250. To convert to the fractional range of 0.0 to 1.0, divide
% the scaled data in the file by 250. To convert to percentage values,
% divide the scaled data in the file by 2.5. Data files may contain
% integers from 0 to 255, as described in the Table 3.
% 
% 0-250   Sea ice concentration (fractional coverage scaled by 250)
% 251 	 Circular mask used in the Arctic to cover the irregularly-shaped
%         data gap around the pole (caused by the orbit inclination and instrument swath) 
% 252     Unused
% 253     Coastlines
% 254     Superimposed land mask
% 255     Missing data 

% -----References:-----
% The plotting code is modified after the documentation for the M_Map
% Toolbox for Matlab (http://www2.ocgy.ubc.ca/~rich/map.html)

% -----Author:-----
% by David L. Jones, July-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2012: now plots imported land mask instead of high resolution
%           coastline; values > 100% are converted to 0% for land masks to
%           these are properly handled by f_NSIDCinterp

% -----Check input & set defaults:-----
% Check that input file is in your path:
if (exist(fname,'file')~=2)
   error('Input file is not in Matlab''s search path!')
end
% -------------------------------------

% Import data:
sfid      = fopen(fname,'r');               % open for reading
txt       = fread(sfid,300,'uint8=>char')'; % read header
P         = fread(sfid,332*316);            % read sea ice values
P         = reshape(P,316,332)';            % reshape into a matrix
P(P==255) = NaN;                            % missing data
P        = (P./250)*100;                    % convert to 0-100 range

% -----Create plot:-----
% From M_Map docummentation: "SSM/I ice products are in a polar
% stereographic projection and the corner points of the grid are given.
% Here we just use those given corner points and 'assume' everything will
% work. It's not bad, although their projection actually uses an
% ellipsoidal earth (m_map uses a spherical earth)."

% M_Map setup:
figure;
set(gcf,'color','w'); % set bg color so printed/exported lakes will be white
m_proj('stereographic','long',0,'lat',-90,'radius',52);

% "Convert bottom and left corner points to screen coords. This
% is of course a kludge" (from M_Map documentation):
[MAPX,nul] = m_ll2xy([225.00 135],[-41.45 -41.45],'clip','off');
[nul,MAPY] = m_ll2xy([317.76 225],[-39.23 -41.45],'clip','off');

% Plot data as an image:
Hi = image(MAPX,MAPY,P);
set(gca,'ydir','normal'); % undo the 'image' command's reversal of ydir

% Customize plot:
colormap([jet(100);0 0 0;1 1 1]);
% m_coast('patch',[.7 .7 .7],'edgecolor','none'); % M_Map's coasline
m_grid('xtick',12,'tickdir','out','ytick',[-80:10:-50], 'XaxisLocation','top');
title({'NSIDC Sea Ice Cover';' '},'fontsize',14,'fontweight','bold');
h = colorbar('v');
set(get(h,'ylabel'),'string','Total Ice Concentration (%)');
   
% Get map coordinates of sea ice data:
Xlim = get(Hi,'Xdata');
Ylim = get(Hi,'Ydata');

% Expand coordinates:
X = linspace(Xlim(1),Xlim(2),size(P,2));
Y = linspace(Ylim(1),Ylim(2),size(P,1));

% After plotting, convert masks to 0% sea ice coverage:
P(P>100) = 0;

% Wrap results up into a structure:
result.txt = txt;
result.P   = P;
result.X   = repmat(X(:)',size(P,1),1);
result.Y   = repmat(Y(:),1,size(P,2));
